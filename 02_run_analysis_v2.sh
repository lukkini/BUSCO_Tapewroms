#!/usr/bin/env bash
# =============================================================================
# 02_run_analysis.sh
# -----------------------------------------------------------------------------
# Core analysis pipeline for building a custom Cestoda BUSCO lineage dataset.
#
# Stages:
#   1. Sequence quality control and cleaning
#   2. Longest-isoform filtering
#   3. FASTA identifier audit (mandatory contract checkpoint)
#   4. OrthoFinder orthogroup inference
#   5. Cestoda marker selection
#   6. Multiple sequence alignment (MAFFT)
#   7. HMM profile construction (hmmbuild)
#   8. Marker-vs-reference hmmsearch tables for downstream score cutoffs
#
# Usage:
#   bash 02_run_analysis.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
shopt -s nullglob

mkdir -p \
    "${CLEAN_PROTEOMES}" \
    "${FILTERED_PROTEOMES}" \
    "${MARKER_DIR}" \
    "${ALIGN_DIR}" \
    "${HMM_DIR}" \
    "${HMMSEARCH_DIR}" \
    "${LOG_DIR}" \
    "${SCRIPTS_DIR}"

LOG="${LOG_DIR}/02_analysis.log"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  02_run_analysis.sh — $(date)"
echo "============================================================"
echo ""

# =============================================================================
# Pre-flight checks
# =============================================================================
N_RAW=$(find "${RAW_PROTEOMES}" -maxdepth 1 -name "*.fa" -size +0c 2>/dev/null | wc -l)
if [[ "${N_RAW}" -eq 0 ]]; then
    echo "[ERROR] No protein FASTA files found in ${RAW_PROTEOMES}/"
    echo "        Run 01_retrieve_data.sh first."
    exit 1
fi

if [[ "${N_RAW}" -lt 26 ]]; then
    echo "[WARN] Expected 26 proteomes but found ${N_RAW}."
    echo "       Continuing, but reduced taxon sampling may affect orthogroup inference."
fi

N_FAIL=$(find "${RAW_PROTEOMES}" -maxdepth 1 -name "*.FAILED" 2>/dev/null | wc -l)
if [[ "${N_FAIL}" -gt 0 ]]; then
    echo "[WARN] ${N_FAIL} species failed to download."
    echo "       Continuing with available proteomes only."
fi

# =============================================================================
# Write helper scripts
# =============================================================================

echo "[INFO] Writing Python helper scripts..."

cat > "${SCRIPTS_DIR}/clean_sequences.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import hashlib
from pathlib import Path
from typing import Dict, List, Set

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def md5(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


def clean_fasta(input_path: str, output_path: str, report_path: str, min_length: int, max_x_frac: float) -> Dict[str, int]:
    seen_hashes: Set[str] = set()
    kept: List[SeqRecord] = []
    counts: Dict[str, int] = {
        "total": 0,
        "too_short": 0,
        "internal_stop": 0,
        "high_x": 0,
        "duplicate": 0,
        "kept": 0,
    }
    report_rows = []

    for rec in SeqIO.parse(input_path, "fasta"):
        counts["total"] += 1
        seq_str = str(rec.seq).upper()
        reason = None

        if len(seq_str) < min_length:
            counts["too_short"] += 1
            reason = "too_short"
        elif "*" in seq_str.rstrip("*"):
            counts["internal_stop"] += 1
            reason = "internal_stop"
        else:
            x_frac = float(seq_str.count("X")) / float(len(seq_str)) if seq_str else 1.0
            if x_frac > max_x_frac:
                counts["high_x"] += 1
                reason = "high_X"

        if reason is None:
            clean_seq = seq_str.rstrip("*")
            h = md5(clean_seq)
            if h in seen_hashes:
                counts["duplicate"] += 1
                reason = "duplicate"
            else:
                seen_hashes.add(h)
                kept.append(SeqRecord(Seq(clean_seq), id=rec.id, name=rec.name, description=rec.description))
                counts["kept"] += 1

        report_rows.append((rec.id, len(seq_str), reason or "kept"))

    with open(output_path, "w") as fh:
        SeqIO.write(kept, fh, "fasta")

    with open(report_path, "w") as fh:
        fh.write("seq_id\tlength\tstatus\n")
        for row in report_rows:
            fh.write("\t".join([str(v) for v in row]) + "\n")

    return counts


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--min_length", type=int, default=30)
    p.add_argument("--max_x_frac", type=float, default=0.10)
    args = p.parse_args()

    counts = clean_fasta(args.input, args.output, args.report, args.min_length, args.max_x_frac)
    pct_kept = 100.0 * counts["kept"] / max(counts["total"], 1)
    print(
        f"  {Path(args.input).name:<50}  "
        f"{counts['total']:>7} in  ->  {counts['kept']:>7} kept  ({pct_kept:.1f}%)  "
        f"[short:{counts['too_short']} stop:{counts['internal_stop']} X:{counts['high_x']} dup:{counts['duplicate']}]"
    )


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/clean_sequences.py"

cat > "${SCRIPTS_DIR}/keep_longest_isoform.py" <<'PYEOF'
#!/usr/bin/env python3
import re
import sys
from typing import Dict, Tuple

from Bio import SeqIO

ISOFORM_PATTERNS = [
    re.compile(r'\.\d+$'),
    re.compile(r'[-_]m?[Rr][Nn][Aa][-_]?\d+$'),
    re.compile(r'[-_][RT][A-Z]$'),
    re.compile(r'[-_][Tt]ranscript[-_]?\d+$'),
    re.compile(r'[-_][Pp]\d+$'),
]


def to_gene_id(seq_id: str) -> str:
    gene_id = seq_id
    for pat in ISOFORM_PATTERNS:
        stripped = pat.sub("", gene_id)
        if stripped != gene_id:
            return stripped
    return gene_id


def keep_longest(input_fasta: str, output_fasta: str) -> Tuple[int, int]:
    genes: Dict[str, object] = {}
    n_in = 0
    for rec in SeqIO.parse(input_fasta, "fasta"):
        n_in += 1
        original_id = rec.id
        gene_id = to_gene_id(original_id)
        if gene_id not in genes or len(rec.seq) > len(genes[gene_id].seq):
            rec.description = f"original_id={original_id} {rec.description}".strip()
            rec.id = gene_id
            rec.name = gene_id
            genes[gene_id] = rec

    with open(output_fasta, "w") as out:
        SeqIO.write(list(genes.values()), out, "fasta")

    return n_in, len(genes)


def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: keep_longest_isoform.py <input.fa> <output.fa>")
        sys.exit(1)

    n_in, n_out = keep_longest(sys.argv[1], sys.argv[2])
    removed = n_in - n_out
    pct = 100.0 * n_out / max(n_in, 1)
    print(
        f"  {sys.argv[1]:<55}  {n_in:>7} isoforms  ->  {n_out:>7} genes  "
        f"({removed} isoforms removed, {pct:.1f}% retained)"
    )


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/keep_longest_isoform.py"

cat > "${SCRIPTS_DIR}/audit_fasta_headers.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Dict, List, Set

from Bio import SeqIO

ILLEGAL_RE = re.compile(r'[\s|:]')
VALID_ID_RE = re.compile(r'^[A-Za-z0-9_.-]+$')


def audit_fasta(path: Path) -> Dict[str, int]:
    seen: Set[str] = set()
    total = 0
    empty = 0
    dup = 0
    illegal = 0
    for rec in SeqIO.parse(str(path), 'fasta'):
        total += 1
        rid = str(rec.id).strip()
        if not rid:
            empty += 1
            continue
        if rid in seen:
            dup += 1
        seen.add(rid)
        if ILLEGAL_RE.search(rid) or not VALID_ID_RE.fullmatch(rid):
            illegal += 1
    return {"total": total, "empty": empty, "dup": dup, "illegal": illegal}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input_dir", required=True)
    p.add_argument("--report", required=True)
    args = p.parse_args()

    input_dir = Path(args.input_dir)
    paths = sorted(input_dir.glob("*.fa"))
    if not paths:
        raise SystemExit(f"[ERROR] No FASTA files found in {input_dir}")

    bad_any = False
    rows: List[List[str]] = []
    for path in paths:
        stats = audit_fasta(path)
        status = "ok" if stats["empty"] == 0 and stats["dup"] == 0 and stats["illegal"] == 0 else "FAIL"
        if status != "ok":
            bad_any = True
        rows.append([path.name, str(stats["total"]), str(stats["empty"]), str(stats["dup"]), str(stats["illegal"]), status])
        print(f"  {path.name:<40} {stats['total']:>7} seqs  empty={stats['empty']}  dup={stats['dup']}  illegal={stats['illegal']}  {status}")

    with open(args.report, "w") as fh:
        fh.write("file\ttotal\tempty_ids\tduplicate_ids\tillegal_ids\tstatus\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    if bad_any:
        raise SystemExit("[ERROR] FASTA identifier audit failed. See report for details.")
    print("\n[OK] FASTA identifier audit passed.")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/audit_fasta_headers.py"

cat > "${SCRIPTS_DIR}/filter_cestoda_markers.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import math
import os
import shutil
from pathlib import Path
from typing import Dict, List

import pandas as pd

CESTODA = frozenset({
    "echinococcus_canadensis_g7",
    "echinococcus_granulosus_g1",
    "echinococcus_multilocularis",
    "hymenolepis_diminuta",
    "hymenolepis_microstoma",
    "hymenolepis_nana",
    "mesocestoides_corti",
    "taenia_asiatica",
    "taenia_multiceps",
    "taenia_saginata",
    "taenia_solium",
})

OUTGROUPS = frozenset({
    "clonorchis_sinensis",
    "fasciola_hepatica",
    "heterobilharzia_americana",
    "opisthorchis_felineus",
    "opisthorchis_viverrini",
    "paragonimus_westermani",
    "schistosoma_mansoni",
    "schistosoma_haematobium",
    "schistosoma_japonicum",
    "schistosoma_rodhaini",
    "trichobilharzia_regenti",
    "trichobilharzia_szidati",
    "schmidtea_mediterranea",
    "gyrodactylus_bullatarudis",
    "gyrodactylus_salaris",
})


def count_copies(cell) -> int:
    if pd.isna(cell):
        return 0
    value = str(cell).strip()
    if value == "" or value.lower() in {"nan", "none"}:
        return 0
    return len([x for x in value.split(",") if x.strip()])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--og_table", required=True)
    p.add_argument("--og_seqdir", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--min_presence", type=float, default=0.80)
    p.add_argument("--min_singlecopy", type=float, default=0.90)
    p.add_argument("--outgroup_min", type=int, default=3)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    df = pd.read_csv(args.og_table, sep="\t", index_col=0)
    cestoda_cols = [c for c in df.columns if c in CESTODA]
    outgroup_cols = [c for c in df.columns if c in OUTGROUPS]
    if not cestoda_cols:
        raise SystemExit("[ERROR] No Cestoda columns found in Orthogroups.tsv")

    n_cestoda = len(cestoda_cols)
    presence_min = int(math.ceil(n_cestoda * args.min_presence))

    kept: List[str] = []
    metadata: List[Dict[str, object]] = []
    removed_presence = 0
    removed_multicopy = 0

    for og_id, row in df.iterrows():
        copies = {sp: count_copies(row[sp]) for sp in cestoda_cols}
        n_present = sum(1 for c in copies.values() if c >= 1)
        n_single_copy = sum(1 for c in copies.values() if c == 1)

        if n_present < presence_min:
            removed_presence += 1
            continue

        singlecopy_fraction = float(n_single_copy) / float(n_present)
        if singlecopy_fraction < args.min_singlecopy:
            removed_multicopy += 1
            continue

        n_outgroup = sum(1 for sp in outgroup_cols if count_copies(row[sp]) >= 1)
        conservation = "platyhelminthes_conserved" if n_outgroup >= args.outgroup_min else "cestoda_specific"
        kept.append(og_id)
        metadata.append({
            "orthogroup_id": og_id,
            "n_cestoda_present": n_present,
            "n_cestoda_singlecopy": n_single_copy,
            "singlecopy_fraction": round(singlecopy_fraction, 4),
            "n_outgroup_present": n_outgroup,
            "conservation": conservation,
        })

    pd.DataFrame(metadata).to_csv(os.path.join(args.out_dir, "marker_metadata.tsv"), sep="\t", index=False)

    seqdir = Path(args.og_seqdir)
    copied = 0
    for og in kept:
        src = seqdir / f"{og}.fa"
        if not src.exists():
            raise SystemExit(f"[ERROR] Missing orthogroup FASTA for retained marker: {src}")
        shutil.copy(src, os.path.join(args.out_dir, f"{og}.fa"))
        copied += 1

    print(f"  Total orthogroups      : {len(df)}")
    print(f"  Removed low presence   : {removed_presence}")
    print(f"  Removed multicopy      : {removed_multicopy}")
    print(f"  Retained markers       : {len(kept)}")
    print(f"  Marker FASTAs copied   : {copied}")

    if copied == 0:
        raise SystemExit("[ERROR] No marker FASTAs were copied")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/filter_cestoda_markers.py"

# =============================================================================
# Stage 1 — Sequence quality control
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Sequence quality control"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

QC_LOG="${LOG_DIR}/qc_summary.tsv"
printf 'species\tn_raw\tn_clean\tpct_kept\tn_short\tn_internal_stop\tn_high_x\tn_duplicate\n' > "${QC_LOG}"

for RAW_FA in "${RAW_PROTEOMES}"/*.fa; do
    TAG="$(basename "${RAW_FA}" .fa)"
    CLEAN_FA="${CLEAN_PROTEOMES}/${TAG}.fa"
    REPORT="${CLEAN_PROTEOMES}/${TAG}.qc.tsv"

    if [[ -s "${CLEAN_FA}" && -s "${REPORT}" ]]; then
        echo "  [SKIP] ${TAG}"
        continue
    fi

    rm -f "${CLEAN_FA}" "${REPORT}"
    python "${SCRIPTS_DIR}/clean_sequences.py" \
        --input "${RAW_FA}" \
        --output "${CLEAN_FA}" \
        --report "${REPORT}" \
        --min_length "${MIN_PROTEIN_LENGTH}" \
        --max_x_frac "${MAX_X_FRACTION}"
done

for REPORT in "${CLEAN_PROTEOMES}"/*.qc.tsv; do
    [[ -f "${REPORT}" ]] || continue
    TAG="$(basename "${REPORT}" .qc.tsv)"
    N_RAW_TAG=$(awk 'NR>1{n++} END{print n+0}' "${REPORT}")
    N_CLEAN_TAG=$(awk -F'\t' 'NR>1 && $3=="kept"{n++} END{print n+0}' "${REPORT}")
    N_SHORT=$(awk -F'\t' 'NR>1 && $3=="too_short"{n++} END{print n+0}' "${REPORT}")
    N_STOP=$(awk -F'\t' 'NR>1 && $3=="internal_stop"{n++} END{print n+0}' "${REPORT}")
    N_X=$(awk -F'\t' 'NR>1 && $3=="high_X"{n++} END{print n+0}' "${REPORT}")
    N_DUP=$(awk -F'\t' 'NR>1 && $3=="duplicate"{n++} END{print n+0}' "${REPORT}")
    PCT=$(awk -v kept="${N_CLEAN_TAG}" -v raw="${N_RAW_TAG}" 'BEGIN{if(raw==0) print "0.00"; else printf "%.2f", (100*kept/raw)}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TAG}" "${N_RAW_TAG}" "${N_CLEAN_TAG}" "${PCT}" "${N_SHORT}" "${N_STOP}" "${N_X}" "${N_DUP}" \
        >> "${QC_LOG}"
done

echo "  QC summary: ${QC_LOG}"

# =============================================================================
# Stage 2 — Longest-isoform filtering
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 2 — Longest-isoform filtering"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

for CLEAN_FA in "${CLEAN_PROTEOMES}"/*.fa; do
    TAG="$(basename "${CLEAN_FA}" .fa)"
    FILT_FA="${FILTERED_PROTEOMES}/${TAG}.fa"

    if [[ -s "${FILT_FA}" ]]; then
        echo "  [SKIP] ${TAG}"
        continue
    fi

    rm -f "${FILT_FA}"
    python "${SCRIPTS_DIR}/keep_longest_isoform.py" "${CLEAN_FA}" "${FILT_FA}"
done

N_FILT=$(find "${FILTERED_PROTEOMES}" -maxdepth 1 -name "*.fa" -size +0c | wc -l)
echo ""
echo "  Filtered proteomes: ${N_FILT} files in ${FILTERED_PROTEOMES}/"

# =============================================================================
# Stage 3 — FASTA identifier audit
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 3 — FASTA identifier audit"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

HEADER_AUDIT_REPORT="${LOG_DIR}/filtered_proteome_header_audit.tsv"
python "${SCRIPTS_DIR}/audit_fasta_headers.py" \
    --input_dir "${FILTERED_PROTEOMES}" \
    --report "${HEADER_AUDIT_REPORT}"

echo "  Header audit report: ${HEADER_AUDIT_REPORT}"

# =============================================================================
# Stage 4 — OrthoFinder
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 4 — OrthoFinder"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

OF_RUN_NAME="cestoda_run1"
OF_OUTPUT_DIR="${OG_DIR}/${OF_RUN_NAME}_output"
OF_RESULTS_DIR="${OF_OUTPUT_DIR}/Results_${OF_RUN_NAME}"

ORTHOGROUP_TABLE_CANDIDATE="$(find "${OF_RESULTS_DIR}" -name "Orthogroups.tsv" -type f 2>/dev/null | head -1 || true)"
ORTHOGROUP_SEQDIR_CANDIDATE="$(find "${OF_RESULTS_DIR}" -name "Orthogroup_Sequences" -type d 2>/dev/null | head -1 || true)"
if [[ -n "${ORTHOGROUP_TABLE_CANDIDATE}" && -f "${ORTHOGROUP_TABLE_CANDIDATE}" && -n "${ORTHOGROUP_SEQDIR_CANDIDATE}" && -d "${ORTHOGROUP_SEQDIR_CANDIDATE}" ]]; then
    echo "  [SKIP] Existing complete OrthoFinder results found: ${OF_RESULTS_DIR}"
else
    if [[ -e "${OF_OUTPUT_DIR}" ]]; then
        echo "  [INFO] Removing stale OrthoFinder output directory: ${OF_OUTPUT_DIR}"
        rm -rf "${OF_OUTPUT_DIR}"
    fi

    orthofinder \
        -f "${FILTERED_PROTEOMES}" \
        -S diamond_ultra_sens \
        -M msa \
        -A mafft \
        -T fasttree \
        -t "${OF_THREADS}" \
        -a "${OF_ANALYSIS_THREADS}" \
        -n "${OF_RUN_NAME}" \
        -o "${OF_OUTPUT_DIR}" \
        2>&1 | tee "${LOG_DIR}/orthofinder.log"
fi

if [[ ! -d "${OF_RESULTS_DIR}" ]]; then
    echo "[ERROR] OrthoFinder results directory not found: ${OF_RESULTS_DIR}"
    exit 1
fi

OG_TABLE="$(find "${OF_RESULTS_DIR}" -name "Orthogroups.tsv" -type f | head -1 || true)"
OG_SEQDIR="$(find "${OF_RESULTS_DIR}" -name "Orthogroup_Sequences" -type d | head -1 || true)"

if [[ -z "${OG_TABLE}" || ! -f "${OG_TABLE}" ]]; then
    echo "[ERROR] Orthogroups.tsv not found under ${OF_RESULTS_DIR}/"
    exit 1
fi
if [[ -z "${OG_SEQDIR}" || ! -d "${OG_SEQDIR}" ]]; then
    echo "[ERROR] Orthogroup_Sequences/ not found under ${OF_RESULTS_DIR}/"
    exit 1
fi

N_ORTHOGROUPS=$(( $(wc -l < "${OG_TABLE}") - 1 ))
echo "  Orthogroup table : ${OG_TABLE}"
echo "  Sequence dir     : ${OG_SEQDIR}"
echo "  Total orthogroups: ${N_ORTHOGROUPS}"

# =============================================================================
# Stage 5 — Marker selection
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 5 — Marker selection"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

MARKER_META="${MARKER_DIR}/marker_metadata.tsv"
EXISTING_MARKERS=$(find "${MARKER_DIR}" -maxdepth 1 -name "OG*.fa" -size +0c | wc -l)

if [[ -s "${MARKER_META}" && "${EXISTING_MARKERS}" -gt 0 ]]; then
    echo "  [SKIP] Existing marker FASTAs found in ${MARKER_DIR}/"
else
    rm -f "${MARKER_DIR}"/OG*.fa "${MARKER_META}"
    python "${SCRIPTS_DIR}/filter_cestoda_markers.py" \
        --og_table "${OG_TABLE}" \
        --og_seqdir "${OG_SEQDIR}" \
        --out_dir "${MARKER_DIR}" \
        --min_presence "${MIN_CESTODA_PRESENCE}" \
        --min_singlecopy "${MIN_SINGLECOPY_FRAC}" \
        --outgroup_min "${OUTGROUP_CONSERVATION_MIN}" \
        2>&1 | tee "${LOG_DIR}/filter_markers.log"
fi

N_MARKERS=$(find "${MARKER_DIR}" -maxdepth 1 -name "OG*.fa" -size +0c | wc -l)
echo "  Markers retained: ${N_MARKERS}"
if [[ "${N_MARKERS}" -eq 0 ]]; then
    echo "[ERROR] No marker FASTAs produced."
    exit 1
fi

# =============================================================================
# Stage 6 — Multiple sequence alignment (MAFFT)
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 6 — Multiple sequence alignment (MAFFT)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Jobs: ${PARALLEL_JOBS} parallel x ${THREADS_PER_JOB} threads each"
echo ""

N_ALIGN_TOTAL=0
run_mafft_pool() {
    local pids=()
    local n_new=0
    local n_skip=0
    local n_fail=0
    local marker_fa og msa

    for marker_fa in "${MARKER_DIR}"/OG*.fa; do
        [[ -f "${marker_fa}" ]] || continue
        og="$(basename "${marker_fa}" .fa)"
        msa="${ALIGN_DIR}/${og}.msa"

        if [[ -s "${msa}" ]]; then
            n_skip=$((n_skip + 1))
            continue
        fi

        rm -f "${msa}"
        mafft --auto --reorder --thread "${THREADS_PER_JOB}" "${marker_fa}" > "${msa}" 2>>"${LOG_DIR}/mafft.log" &
        pids+=("$!")
        n_new=$((n_new + 1))

        if [[ "${#pids[@]}" -ge "${PARALLEL_JOBS}" ]]; then
            for pid in "${pids[@]}"; do
                if ! wait "${pid}"; then
                    n_fail=$((n_fail + 1))
                fi
            done
            pids=()
            echo "  Aligned ${n_new} marker families so far..."
        fi
    done

    for pid in "${pids[@]}"; do
        if ! wait "${pid}"; then
            n_fail=$((n_fail + 1))
        fi
    done

    N_ALIGN_TOTAL=$(find "${ALIGN_DIR}" -maxdepth 1 -name "*.msa" -size +0c | wc -l)
    echo "  Alignments: ${n_new} new / ${n_skip} existing / ${n_fail} failed"
    echo "  Total MSA files: ${N_ALIGN_TOTAL}"

    if [[ "${n_fail}" -ne 0 ]]; then
        echo "[ERROR] One or more MAFFT jobs failed. See ${LOG_DIR}/mafft.log"
        exit 1
    fi
}

run_mafft_pool
if [[ "${N_ALIGN_TOTAL}" -ne "${N_MARKERS}" ]]; then
    echo "[ERROR] Alignment count (${N_ALIGN_TOTAL}) does not match marker count (${N_MARKERS})."
    exit 1
fi

# =============================================================================
# Stage 7 — HMM profile construction
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 7 — HMM profile construction (hmmbuild)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Jobs: ${PARALLEL_JOBS} parallel x ${THREADS_PER_JOB} threads each"
echo ""

run_hmmbuild_pool() {
    local pids=()
    local n_new=0
    local n_skip=0
    local n_fail=0
    local msa og hmm

    for msa in "${ALIGN_DIR}"/*.msa; do
        [[ -f "${msa}" ]] || continue
        og="$(basename "${msa}" .msa)"
        hmm="${HMM_DIR}/${og}.hmm"

        if [[ -s "${hmm}" ]]; then
            n_skip=$((n_skip + 1))
            continue
        fi

        rm -f "${hmm}"
        hmmbuild --cpu "${THREADS_PER_JOB}" --amino -n "${og}" "${hmm}" "${msa}" >> "${LOG_DIR}/hmmbuild.log" 2>&1 &
        pids+=("$!")
        n_new=$((n_new + 1))

        if [[ "${#pids[@]}" -ge "${PARALLEL_JOBS}" ]]; then
            for pid in "${pids[@]}"; do
                if ! wait "${pid}"; then
                    n_fail=$((n_fail + 1))
                fi
            done
            pids=()
            echo "  Built ${n_new} HMMs so far..."
        fi
    done

    for pid in "${pids[@]}"; do
        if ! wait "${pid}"; then
            n_fail=$((n_fail + 1))
        fi
    done

    echo "  hmmbuild: ${n_new} new / ${n_skip} existing / ${n_fail} failed"
    if [[ "${n_fail}" -ne 0 ]]; then
        echo "[ERROR] One or more hmmbuild jobs failed. See ${LOG_DIR}/hmmbuild.log"
        exit 1
    fi
}

run_hmmbuild_pool
N_HMM_TOTAL=$(find "${HMM_DIR}" -maxdepth 1 -name "*.hmm" -size +0c | wc -l)
echo "  Total HMM profiles: ${N_HMM_TOTAL}"
if [[ "${N_HMM_TOTAL}" -ne "${N_MARKERS}" ]]; then
    echo "[ERROR] HMM count (${N_HMM_TOTAL}) does not match marker count (${N_MARKERS})."
    exit 1
fi

# =============================================================================
# Stage 8 — Marker-vs-reference hmmsearch tables
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 8 — hmmsearch reference tables"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Jobs: ${PARALLEL_JOBS} parallel x ${THREADS_PER_JOB} threads each"
echo ""

ALL_PROTEINS="${HMMSEARCH_DIR}/all_reference_proteins.faa"
if [[ ! -s "${ALL_PROTEINS}" ]]; then
    cat "${FILTERED_PROTEOMES}"/*.fa > "${ALL_PROTEINS}"
fi
if [[ ! -s "${ALL_PROTEINS}" ]]; then
    echo "[ERROR] Failed to build concatenated reference protein database: ${ALL_PROTEINS}"
    exit 1
fi

run_hmmsearch_pool() {
    local pids=()
    local n_new=0
    local n_skip=0
    local n_fail=0
    local hmm og tbl out

    for hmm in "${HMM_DIR}"/*.hmm; do
        [[ -f "${hmm}" ]] || continue
        og="$(basename "${hmm}" .hmm)"
        tbl="${HMMSEARCH_DIR}/${og}.tbl"
        out="${HMMSEARCH_DIR}/${og}.hmmsearch.out"

        if [[ -s "${tbl}" && "${tbl}" -nt "${hmm}" && "${tbl}" -nt "${ALL_PROTEINS}" ]]; then
            n_skip=$((n_skip + 1))
            continue
        fi

        rm -f "${tbl}" "${out}"
        hmmsearch --cpu "${THREADS_PER_JOB}" --noali --tblout "${tbl}" "${hmm}" "${ALL_PROTEINS}" > "${out}" 2>>"${LOG_DIR}/hmmsearch.log" &
        pids+=("$!")
        n_new=$((n_new + 1))

        if [[ "${#pids[@]}" -ge "${PARALLEL_JOBS}" ]]; then
            for pid in "${pids[@]}"; do
                if ! wait "${pid}"; then
                    n_fail=$((n_fail + 1))
                fi
            done
            pids=()
            echo "  Generated ${n_new} hmmsearch tables so far..."
        fi
    done

    for pid in "${pids[@]}"; do
        if ! wait "${pid}"; then
            n_fail=$((n_fail + 1))
        fi
    done

    echo "  hmmsearch tables: ${n_new} new / ${n_skip} existing / ${n_fail} failed"
    if [[ "${n_fail}" -ne 0 ]]; then
        echo "[ERROR] One or more hmmsearch jobs failed. See ${LOG_DIR}/hmmsearch.log"
        exit 1
    fi
}

run_hmmsearch_pool
N_TBL_TOTAL=$(find "${HMMSEARCH_DIR}" -maxdepth 1 -name "*.tbl" -size +0c | wc -l)
echo "  Total hmmsearch tables: ${N_TBL_TOTAL}"
if [[ "${N_TBL_TOTAL}" -ne "${N_HMM_TOTAL}" ]]; then
    echo "[ERROR] hmmsearch table count (${N_TBL_TOTAL}) does not match HMM count (${N_HMM_TOTAL})."
    exit 1
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Analysis pipeline complete — $(date)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Species processed        : ${N_FILT}"
echo "  Orthogroups inferred     : ${N_ORTHOGROUPS}"
echo "  Markers selected         : ${N_MARKERS}"
echo "  Alignments produced      : ${N_ALIGN_TOTAL}"
echo "  HMM profiles built       : ${N_HMM_TOTAL}"
echo "  hmmsearch tables built   : ${N_TBL_TOTAL}"
echo ""
echo "  Log: ${LOG}"
echo "  Next step: bash 03_build_odb.sh"
echo ""
