#!/usr/bin/env bash
# =============================================================================
# 03_build_odb.sh
# -----------------------------------------------------------------------------
# Export a functionally valid custom BUSCO v6 lineage dataset for Cestoda.
#
# Critical contract rule:
#   The FINAL exported lineage uses ONE numeric BUSCO marker ID namespace.
#   That final numeric ID must match consistently across:
#     - HMM filenames
#     - internal HMM NAME fields
#     - scores_cutoff IDs
#     - info/ogs.id.info IDs
#     - links_to_ODB12.txt IDs
#     - refseq_db.faa.gz FASTA primary IDs
#     - ancestral / ancestral_variants headers
#
# Usage:
#   bash 03_build_odb.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

if [[ "${DATASET_NAME}" == *_odb12 ]]; then
    FINAL_DATASET_NAME="${DATASET_NAME}"
else
    FINAL_DATASET_NAME="${DATASET_NAME}_odb12"
fi
FINAL_ODB_DIR="${PROJECT_ROOT}/${FINAL_DATASET_NAME}"

mkdir -p "${LOG_DIR}" "${SCRIPTS_DIR}" "${VALIDATION_DIR}"

LOG="${LOG_DIR}/03_build_odb.log"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  03_build_odb.sh — $(date)"
echo "============================================================"
echo ""
echo "[INFO] Final lineage name : ${FINAL_DATASET_NAME}"
echo "[INFO] Final lineage dir  : ${FINAL_ODB_DIR}"
echo ""

# =============================================================================
# Runtime resources
# =============================================================================
if command -v nproc >/dev/null 2>&1; then
    AVAILABLE_CORES="$(nproc)"
else
    AVAILABLE_CORES="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi

if [[ -z "${AVAILABLE_CORES}" ]] || ! [[ "${AVAILABLE_CORES}" =~ ^[0-9]+$ ]] || [[ "${AVAILABLE_CORES}" -lt 1 ]]; then
    AVAILABLE_CORES=1
fi

echo "[INFO] Available CPU cores : ${AVAILABLE_CORES}"

# =============================================================================
# Pre-flight checks
# =============================================================================
N_HMM=$(find "${HMM_DIR}" -maxdepth 1 -name "*.hmm" -size +0c 2>/dev/null | wc -l)
if [[ "${N_HMM}" -eq 0 ]]; then
    echo "[ERROR] No HMM profiles found in ${HMM_DIR}/"
    echo "        Run 02_run_analysis.sh first."
    exit 1
fi

N_TBL=$(find "${HMMSEARCH_DIR}" -maxdepth 1 -name "*.tbl" -size +0c 2>/dev/null | wc -l)
if [[ "${N_TBL}" -eq 0 ]]; then
    echo "[ERROR] No hmmsearch tables found in ${HMMSEARCH_DIR}/"
    echo "        02_run_analysis.sh must finish Stage 8 successfully first."
    exit 1
fi

N_MARKER_FASTA=$(find "${MARKER_DIR}" -maxdepth 1 -name "OG*.fa" -size +0c 2>/dev/null | wc -l)
if [[ "${N_MARKER_FASTA}" -eq 0 ]]; then
    echo "[ERROR] No marker FASTAs found in ${MARKER_DIR}/"
    exit 1
fi

echo "[INFO] Source HMMs          : ${N_HMM}"
echo "[INFO] Source hmmsearch tbl : ${N_TBL}"
echo "[INFO] Marker FASTAs        : ${N_MARKER_FASTA}"

# =============================================================================
# Fresh export directory
# =============================================================================

echo ""
echo "[INFO] Resetting export-derived files to avoid stale broken state..."
rm -rf "${FINAL_ODB_DIR}"
mkdir -p \
    "${FINAL_ODB_DIR}/hmms" \
    "${FINAL_ODB_DIR}/info" \
    "${FINAL_ODB_DIR}/prfl"

# =============================================================================
# Write helper scripts
# =============================================================================

cat > "${SCRIPTS_DIR}/audit_export_inputs.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Set


def stems_from_dir(path: Path, pattern: str) -> Set[str]:
    return {p.stem for p in path.glob(pattern) if p.is_file() and p.stat().st_size > 0}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--hmm_dir", required=True)
    p.add_argument("--tbl_dir", required=True)
    p.add_argument("--marker_dir", required=True)
    args = p.parse_args()

    hmm_ids = stems_from_dir(Path(args.hmm_dir), "*.hmm")
    tbl_ids = stems_from_dir(Path(args.tbl_dir), "*.tbl")
    marker_ids = stems_from_dir(Path(args.marker_dir), "OG*.fa")

    if not hmm_ids:
        raise SystemExit("[ERROR] No HMMs found for export audit")
    if not tbl_ids:
        raise SystemExit("[ERROR] No hmmsearch tables found for export audit")
    if not marker_ids:
        raise SystemExit("[ERROR] No marker FASTAs found for export audit")

    if hmm_ids != tbl_ids or hmm_ids != marker_ids:
        missing_in_tbl = sorted(hmm_ids - tbl_ids)
        missing_in_marker = sorted(hmm_ids - marker_ids)
        extra_tbl = sorted(tbl_ids - hmm_ids)
        extra_marker = sorted(marker_ids - hmm_ids)

        msg = (
            "[ERROR] Export input contract mismatch across HMMs / tbl / marker FASTAs.\n"
            f"        HMMs={len(hmm_ids)} tbl={len(tbl_ids)} markers={len(marker_ids)}\n"
            f"        Missing tbl for HMM IDs: {missing_in_tbl[:10]}\n"
            f"        Missing marker FASTAs for HMM IDs: {missing_in_marker[:10]}\n"
            f"        Extra tbl IDs without HMM: {extra_tbl[:10]}\n"
            f"        Extra marker IDs without HMM: {extra_marker[:10]}"
        )
        raise SystemExit(msg)

    print(f"  Export input audit passed: {len(hmm_ids)} OG IDs shared across HMMs, tbls, and marker FASTAs.")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/audit_export_inputs.py"

cat > "${SCRIPTS_DIR}/build_busco_id_map.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--hmm_dir", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--prefix", type=int, default=1000000000)
    args = p.parse_args()

    hmm_dir = Path(args.hmm_dir)
    og_ids = sorted([p.stem for p in hmm_dir.glob("*.hmm")])
    if not og_ids:
        raise SystemExit(f"[ERROR] No HMMs found in {hmm_dir}")

    seen = set()
    for og in og_ids:
        if og in seen:
            raise SystemExit(f"[ERROR] Duplicate source HMM stem detected: {og}")
        seen.add(og)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["orthogroup_id", "busco_id"])
        for idx, og in enumerate(og_ids, start=1):
            writer.writerow([og, str(args.prefix + idx)])

    print(f"  Written BUSCO ID map: {args.output} ({len(og_ids)} entries)")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/build_busco_id_map.py"

cat > "${SCRIPTS_DIR}/rewrite_hmm_ids.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def rewrite_hmm(src: Path, dst: Path, busco_id: str) -> None:
    lines = src.read_text().splitlines()
    out_lines = []
    saw_name = False
    for line in lines:
        if line.startswith("NAME"):
            out_lines.append(f"NAME  {busco_id}")
            saw_name = True
        elif line.startswith("ACC"):
            out_lines.append(f"ACC   {busco_id}")
        else:
            out_lines.append(line)
    if not saw_name:
        raise SystemExit(f"[ERROR] HMM missing NAME field: {src}")
    dst.write_text("\n".join(out_lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--id_map", required=True)
    p.add_argument("--src_hmm_dir", required=True)
    p.add_argument("--dst_hmm_dir", required=True)
    args = p.parse_args()

    src_dir = Path(args.src_hmm_dir)
    dst_dir = Path(args.dst_hmm_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    with open(args.id_map, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            og_id = row["orthogroup_id"]
            busco_id = row["busco_id"]
            src = src_dir / f"{og_id}.hmm"
            dst = dst_dir / f"{busco_id}.hmm"
            if not src.exists():
                raise SystemExit(f"[ERROR] Missing source HMM: {src}")
            rewrite_hmm(src, dst, busco_id)
            n += 1

    print(f"  Rewritten/exported HMMs: {n}")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/rewrite_hmm_ids.py"

cat > "${SCRIPTS_DIR}/compute_score_cutoffs_from_map.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from typing import List


def parse_tblout(path: Path) -> List[float]:
    scores: List[float] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            try:
                scores.append(float(cols[5]))
            except ValueError:
                continue
    return scores


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--id_map", required=True)
    p.add_argument("--tbl_dir", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--fraction", type=float, default=0.90)
    args = p.parse_args()

    rows = []
    with open(args.id_map, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append((row["orthogroup_id"], row["busco_id"]))

    if not rows:
        raise SystemExit("[ERROR] Empty BUSCO ID map")

    out_lines = []
    for og_id, busco_id in rows:
        tbl = Path(args.tbl_dir) / f"{og_id}.tbl"
        if not tbl.exists():
            raise SystemExit(f"[ERROR] Missing hmmsearch table: {tbl}")
        scores = parse_tblout(tbl)
        if not scores:
            raise SystemExit(f"[ERROR] Empty/scoreless hmmsearch table for {og_id}: {tbl}")
        cutoff = round(min(scores) * args.fraction, 2)
        out_lines.append((int(busco_id), cutoff))

    out_lines.sort(key=lambda x: x[0])
    with open(args.output, "w") as out:
        for busco_id, cutoff in out_lines:
            out.write(f"{busco_id}\t{cutoff}\n")

    print(f"  Written score cutoffs: {len(out_lines)}")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/compute_score_cutoffs_from_map.py"

cat > "${SCRIPTS_DIR}/build_links_refseq_and_ancestral.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def parse_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    name = None
    seq_parts: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts)
                name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if name is not None:
            yield name, "".join(seq_parts)


def load_metadata(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    out: Dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        og_field = "orthogroup_id" if "orthogroup_id" in fieldnames else None
        desc_field = "conservation" if "conservation" in fieldnames else None
        if og_field is None:
            return out
        for row in reader:
            og = str(row.get(og_field, "")).strip()
            if not og:
                continue
            out[og] = str(row.get(desc_field, "custom_cestoda_marker")).strip() or "custom_cestoda_marker"
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--id_map", required=True)
    p.add_argument("--marker_dir", required=True)
    p.add_argument("--meta", required=True)
    p.add_argument("--links_out", required=True)
    p.add_argument("--refseq_out", required=True)
    p.add_argument("--provenance_out", required=True)
    p.add_argument("--ancestral_out", required=True)
    p.add_argument("--ancestral_variants_out", required=True)
    p.add_argument("--src_hmm_dir", required=True)
    args = p.parse_args()

    metadata = load_metadata(Path(args.meta))
    marker_dir = Path(args.marker_dir)
    src_hmm_dir = Path(args.src_hmm_dir)

    rows = []
    with open(args.id_map, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append((row["orthogroup_id"], row["busco_id"]))

    with open(args.links_out, "w") as links_fh, \
         gzip.open(args.refseq_out, "wt") as ref_fh, \
         open(args.provenance_out, "w") as prov_fh, \
         open(args.ancestral_out, "w") as anc_fh, \
         open(args.ancestral_variants_out, "w") as ancv_fh:

        prov_fh.write("busco_id\torthogroup_id\trepresentative_header\trepresentative_length\n")

        for og_id, busco_id in rows:
            marker_fa = marker_dir / f"{og_id}.fa"
            if not marker_fa.exists():
                raise SystemExit(f"[ERROR] Missing marker FASTA: {marker_fa}")

            desc = metadata.get(og_id, "custom_cestoda_marker")
            links_fh.write(f"{busco_id}\t{desc}\thttps://www.orthodb.org/?query={og_id}\n")

            best_header = None
            best_seq = ""
            variant_idx = 0
            for header, seq in parse_fasta(marker_fa):
                if not seq:
                    continue
                variant_idx += 1
                ancv_fh.write(f">{busco_id}_{variant_idx}\n")
                for i in range(0, len(seq), 80):
                    ancv_fh.write(seq[i:i+80] + "\n")
                if len(seq) > len(best_seq):
                    best_header = header
                    best_seq = seq

            if not best_seq:
                raise SystemExit(f"[ERROR] No usable sequences in {marker_fa}")

            anc_fh.write(f">{busco_id}\n")
            for i in range(0, len(best_seq), 80):
                anc_fh.write(best_seq[i:i+80] + "\n")

            clean_header = (best_header or "representative").replace("\t", " ").replace("\n", " ")
            ref_fh.write(f">{busco_id} orthogroup={og_id} representative={clean_header}\n")
            for i in range(0, len(best_seq), 80):
                ref_fh.write(best_seq[i:i+80] + "\n")

            prov_fh.write(f"{busco_id}\t{og_id}\t{clean_header}\t{len(best_seq)}\n")

    print(f"  Written links      : {args.links_out}")
    print(f"  Written refseq DB  : {args.refseq_out}")
    print(f"  Written provenance : {args.provenance_out}")
    print(f"  Written ancestral  : {args.ancestral_out}")
    print(f"  Written variants   : {args.ancestral_variants_out}")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/build_links_refseq_and_ancestral.py"

cat > "${SCRIPTS_DIR}/audit_lineage_contract.py" <<'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from pathlib import Path
from typing import Set


NUMERIC_RE = re.compile(r'^[0-9]+$')


def read_id_set_from_lines(path: Path, first_column_only: bool = False) -> Set[str]:
    ids: Set[str] = set()
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if first_column_only:
                token = line.split("\t")[0]
            else:
                token = line
            ids.add(token)
    return ids


def read_hmm_name(path: Path) -> str:
    with path.open() as fh:
        for line in fh:
            if line.startswith("NAME"):
                return line.split(None, 1)[1].strip()
    raise RuntimeError(f"No NAME field found in {path}")


def read_refseq_ids(path: Path) -> Set[str]:
    ids: Set[str] = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset_dir", required=True)
    args = p.parse_args()

    ds = Path(args.dataset_dir)
    hmms_dir = ds / "hmms"
    scores = ds / "scores_cutoff"
    ogs = ds / "info" / "ogs.id.info"
    links = ds / "links_to_ODB12.txt"
    refseq = ds / "refseq_db.faa.gz"
    id_map = ds / "info" / "busco_id_map.tsv"

    if not hmms_dir.exists():
        raise SystemExit("[ERROR] hmms/ directory missing")

    hmm_ids = set()
    for hmm in sorted(hmms_dir.glob("*.hmm")):
        stem = hmm.stem
        if not NUMERIC_RE.fullmatch(stem):
            raise SystemExit(f"[ERROR] Non-numeric HMM filename stem: {stem}")
        internal_name = read_hmm_name(hmm)
        if internal_name != stem:
            raise SystemExit(f"[ERROR] HMM NAME mismatch: file={stem} internal={internal_name}")
        hmm_ids.add(stem)

    if not hmm_ids:
        raise SystemExit("[ERROR] No final HMMs found")

    score_ids = read_id_set_from_lines(scores, first_column_only=True)
    ogs_ids = read_id_set_from_lines(ogs)
    links_ids = read_id_set_from_lines(links, first_column_only=True)
    refseq_ids = read_refseq_ids(refseq)

    if not (hmm_ids == score_ids == ogs_ids == links_ids == refseq_ids):
        raise SystemExit(
            "[ERROR] Final lineage ID contract mismatch across exported files.\n"
            f"        hmms={len(hmm_ids)} scores={len(score_ids)} ogs={len(ogs_ids)} links={len(links_ids)} refseq={len(refseq_ids)}"
        )

    with open(id_map, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        mapped_busco_ids = {row["busco_id"] for row in reader}
    if mapped_busco_ids != hmm_ids:
        raise SystemExit("[ERROR] busco_id_map.tsv does not match exported numeric BUSCO IDs")

    print(f"  Contract audit passed for {len(hmm_ids)} BUSCO markers.")


if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/audit_lineage_contract.py"

# =============================================================================
# Stage 1 — BUSCO ID map
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Export input audit + BUSCO ID map"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/audit_export_inputs.py" \
    --hmm_dir "${HMM_DIR}" \
    --tbl_dir "${HMMSEARCH_DIR}" \
    --marker_dir "${MARKER_DIR}"

BUSCO_ID_MAP="${FINAL_ODB_DIR}/info/busco_id_map.tsv"
python "${SCRIPTS_DIR}/build_busco_id_map.py" \
    --hmm_dir "${HMM_DIR}" \
    --output "${BUSCO_ID_MAP}"

N_MAP=$(($(wc -l < "${BUSCO_ID_MAP}") - 1))
if [[ "${N_MAP}" -ne "${N_HMM}" ]]; then
    echo "[ERROR] BUSCO ID map entry count (${N_MAP}) does not match source HMM count (${N_HMM})"
    exit 1
fi

# =============================================================================
# Stage 2 — scores_cutoff from OG-named hmmsearch tables
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 2 — Score cutoff computation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/compute_score_cutoffs_from_map.py" \
    --id_map "${BUSCO_ID_MAP}" \
    --tbl_dir "${HMMSEARCH_DIR}" \
    --output "${FINAL_ODB_DIR}/scores_cutoff" \
    --fraction "${SCORE_CUTOFF_FRACTION}"

N_CUTOFFS=$(wc -l < "${FINAL_ODB_DIR}/scores_cutoff")
if [[ "${N_CUTOFFS}" -ne "${N_HMM}" ]]; then
    echo "[ERROR] scores_cutoff entry count (${N_CUTOFFS}) does not match HMM count (${N_HMM})"
    exit 1
fi

# =============================================================================
# Stage 3 — Final HMM export with internal NAME rewrite
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 3 — Final HMM export"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/rewrite_hmm_ids.py" \
    --id_map "${BUSCO_ID_MAP}" \
    --src_hmm_dir "${HMM_DIR}" \
    --dst_hmm_dir "${FINAL_ODB_DIR}/hmms"

N_FINAL_HMMS=$(find "${FINAL_ODB_DIR}/hmms" -maxdepth 1 -name "*.hmm" -size +0c | wc -l)
if [[ "${N_FINAL_HMMS}" -ne "${N_HMM}" ]]; then
    echo "[ERROR] Final exported HMM count (${N_FINAL_HMMS}) does not match source HMM count (${N_HMM})"
    exit 1
fi

# =============================================================================
# Stage 4 — links_to_ODB12.txt, refseq_db.faa.gz, ancestral, variants
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 4 — Lineage sequence assets"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/build_links_refseq_and_ancestral.py" \
    --id_map "${BUSCO_ID_MAP}" \
    --marker_dir "${MARKER_DIR}" \
    --meta "${MARKER_DIR}/marker_metadata.tsv" \
    --links_out "${FINAL_ODB_DIR}/links_to_ODB12.txt" \
    --refseq_out "${FINAL_ODB_DIR}/refseq_db.faa.gz" \
    --provenance_out "${FINAL_ODB_DIR}/info/refseq_provenance.tsv" \
    --ancestral_out "${FINAL_ODB_DIR}/ancestral" \
    --ancestral_variants_out "${FINAL_ODB_DIR}/ancestral_variants" \
    --src_hmm_dir "${HMM_DIR}"

# =============================================================================
# Stage 5 — prfl placeholders and metadata
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 5 — Dataset metadata"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

while IFS=$'\t' read -r OG_ID BUSCO_ID; do
    [[ "${OG_ID}" == "orthogroup_id" ]] && continue
    : > "${FINAL_ODB_DIR}/prfl/${BUSCO_ID}.prfl"
done < "${BUSCO_ID_MAP}"

awk -F'\t' 'NR>1{print $2}' "${BUSCO_ID_MAP}" | sort -n > "${FINAL_ODB_DIR}/info/ogs.id.info"

TODAY=$(date +%Y-%m-%d)
cat > "${FINAL_ODB_DIR}/dataset.cfg" <<CFGEOF
name=${FINAL_DATASET_NAME}
domain=${DATASET_DOMAIN}
creation_date=${TODAY}
number_of_BUSCOs=${N_HMM}
number_of_species=11
CFGEOF

cat > "${FINAL_ODB_DIR}/info/species.info" <<'EOF2'
echinococcus_canadensis_g7	Cestoda	PRJEB8992	target
echinococcus_granulosus_g1	Cestoda	PRJEB121	target
echinococcus_multilocularis	Cestoda	PRJEB122	target
hymenolepis_diminuta	Cestoda	PRJEB30942	target
hymenolepis_microstoma	Cestoda	PRJEB124	target
hymenolepis_nana	Cestoda	PRJEB508	target
mesocestoides_corti	Cestoda	PRJEB510	target
taenia_asiatica	Cestoda	PRJEB532	target
taenia_multiceps	Cestoda	PRJNA307624	target
taenia_saginata	Cestoda	PRJNA71493	target
taenia_solium	Cestoda	PRJNA170813	target
clonorchis_sinensis	Trematoda	PRJNA386618	outgroup
fasciola_hepatica	Trematoda	PRJEB58756	outgroup
heterobilharzia_americana	Trematoda	TD1_PRJEB44434	outgroup
opisthorchis_felineus	Trematoda	PRJNA413383	outgroup
opisthorchis_viverrini	Trematoda	PRJNA222628	outgroup
paragonimus_westermani	Trematoda	PRJNA219632	outgroup
schistosoma_mansoni	Trematoda	PRJEA36577	outgroup
schistosoma_haematobium	Trematoda	PRJNA78265	outgroup
schistosoma_japonicum	Trematoda	PRJNA520774	outgroup
schistosoma_rodhaini	Trematoda	TD1_PRJEB44434	outgroup
trichobilharzia_regenti	Trematoda	PRJEB44434	outgroup
trichobilharzia_szidati	Trematoda	PRJEB44434	outgroup
schmidtea_mediterranea	Rhabditophora	S2F19H1_PRJNA885486	outgroup
gyrodactylus_bullatarudis	Monogenea	PRJNA532341	outgroup
gyrodactylus_salaris	Monogenea	PRJNA244375	outgroup
EOF2

# =============================================================================
# Stage 6 — Contract audit
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 6 — Final lineage contract audit"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/audit_lineage_contract.py" --dataset_dir "${FINAL_ODB_DIR}"

REQUIRED=(
    "hmms"
    "ancestral"
    "ancestral_variants"
    "scores_cutoff"
    "refseq_db.faa.gz"
    "dataset.cfg"
    "info/ogs.id.info"
    "info/species.info"
    "info/busco_id_map.tsv"
    "info/refseq_provenance.tsv"
    "links_to_ODB12.txt"
    "prfl"
)

for item in "${REQUIRED[@]}"; do
    full="${FINAL_ODB_DIR}/${item}"
    if [[ ! -e "${full}" ]]; then
        echo "[ERROR] Missing required export item: ${full}"
        exit 1
    fi
    if [[ ! -d "${full}" && ! -s "${full}" ]]; then
        echo "[ERROR] Empty required export item: ${full}"
        exit 1
    fi
done

echo ""
echo "════════════════════════════════════════════════════════"
echo "  Custom Cestoda BUSCO lineage export complete"
echo "════════════════════════════════════════════════════════"
echo ""
echo "  Dataset location : ${FINAL_ODB_DIR}/"
echo "  BUSCO markers    : ${N_HMM}"
echo "  ID map           : ${BUSCO_ID_MAP}"
echo ""
echo "  Next step: bash 04_test_cestoda_lineage.sh"
echo ""
