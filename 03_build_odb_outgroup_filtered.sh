#!/usr/bin/env bash
# =============================================================================
# 03_build_odb.sh
# -----------------------------------------------------------------------------
# Build a BUSCO v6 official-compatible custom lineage dataset for Cestoda.
#
# This version adds an explicit OUTGROUP PROTEIN FILTER before exporting the
# final lineage. The filter operates on orthogroup membership statistics from
# the original OrthoFinder table, so the final exported marker set is defined by:
#
#   1. Existing selected marker families from 02_run_analysis.sh
#   2. HMM/tbl/marker family contract consistency
#   3. Additional outgroup-aware filtering prior to final export
#
# Outgroups are therefore used as actual outgroups at export time:
#   - ingroup (Cestoda) defines the lineage core
#   - outgroups are used as a negative / contrast filter
#
# Heavy upstream jobs are NOT recomputed here:
#   - no OrthoFinder
#   - no MAFFT
#   - no hmmbuild
#   - no hmmsearch
#
# This script only reuses their outputs.
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

mkdir -p     "${HMMSEARCH_DIR}"     "${FINAL_ODB_DIR}/hmms"     "${FINAL_ODB_DIR}/info"     "${FINAL_ODB_DIR}/prfl"     "${VALIDATION_DIR}"     "${LOG_DIR}"     "${SCRIPTS_DIR}"

LOG="${LOG_DIR}/03_build_odb.log"
exec > >(tee -a "${LOG}") 2>&1

if command -v nproc &>/dev/null; then
    AVAILABLE_CORES="$(nproc)"
else
    AVAILABLE_CORES="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi

if [[ -z "${AVAILABLE_CORES}" ]] || ! [[ "${AVAILABLE_CORES}" =~ ^[0-9]+$ ]] || [[ "${AVAILABLE_CORES}" -lt 1 ]]; then
    AVAILABLE_CORES=1
fi

echo ""
echo "============================================================"
echo "  03_build_odb.sh — $(date)"
echo "============================================================"
echo ""
echo "[INFO] Final lineage name : ${FINAL_DATASET_NAME}"
echo "[INFO] Final lineage dir  : ${FINAL_ODB_DIR}"
echo ""
echo "[INFO] Available CPU cores : ${AVAILABLE_CORES}"

guess_orthofinder_results_dir() {
    local CAND
    for CAND in         "${OG_DIR}/cestoda_run1_output/Results_cestoda_run1"         "${OG_DIR}/Results_cestoda_run1"
    do
        if [[ -d "${CAND}" ]] && [[ -s "${CAND}/Orthogroups/Orthogroups.tsv" ]]; then
            printf '%s\n' "${CAND}"
            return 0
        fi
    done
    CAND="$(find "${OG_DIR}" -maxdepth 3 -type f -path '*/Orthogroups/Orthogroups.tsv' | sort | head -n 1 || true)"
    if [[ -n "${CAND}" ]]; then
        dirname "$(dirname "${CAND}")"
        return 0
    fi
    return 1
}

ORTHOFINDER_RESULTS_DIR="$(guess_orthofinder_results_dir || true)"
if [[ -z "${ORTHOFINDER_RESULTS_DIR}" ]] || [[ ! -d "${ORTHOFINDER_RESULTS_DIR}" ]]; then
    echo "[ERROR] Could not locate OrthoFinder results directory under ${OG_DIR}/"
    echo "        Run 02_run_analysis.sh first."
    exit 1
fi

ORTHOGROUPS_TSV="${ORTHOFINDER_RESULTS_DIR}/Orthogroups/Orthogroups.tsv"
if [[ ! -s "${ORTHOGROUPS_TSV}" ]]; then
    echo "[ERROR] Missing Orthogroups.tsv: ${ORTHOGROUPS_TSV}"
    exit 1
fi

N_HMM=$(find "${HMM_DIR}" -maxdepth 1 -name "*.hmm" -size +0c 2>/dev/null | wc -l)
N_TBL=$(find "${HMMSEARCH_DIR}" -maxdepth 1 -name "*.tbl" -size +0c 2>/dev/null | wc -l)
N_MARKER_FASTA=$(find "${MARKER_DIR}" -maxdepth 1 -name "OG*.fa" -size +0c 2>/dev/null | wc -l)

if [[ "${N_HMM}" -eq 0 ]]; then
    echo "[ERROR] No HMM profiles found in ${HMM_DIR}/"
    exit 1
fi
if [[ "${N_TBL}" -eq 0 ]]; then
    echo "[ERROR] No hmmsearch tables found in ${HMMSEARCH_DIR}/"
    exit 1
fi
if [[ "${N_MARKER_FASTA}" -eq 0 ]]; then
    echo "[ERROR] No marker FASTA files found in ${MARKER_DIR}/"
    exit 1
fi

echo "[INFO] Source HMMs          : ${N_HMM}"
echo "[INFO] Source hmmsearch tbl : ${N_TBL}"
echo "[INFO] Marker FASTAs        : ${N_MARKER_FASTA}"
echo ""

echo "[INFO] Resetting export-derived files to avoid stale broken state..."
rm -rf "${FINAL_ODB_DIR}"
mkdir -p "${FINAL_ODB_DIR}/hmms" "${FINAL_ODB_DIR}/info" "${FINAL_ODB_DIR}/prfl"

cat > "${SCRIPTS_DIR}/audit_export_inputs.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
from pathlib import Path

def stems_from_dir(path_str: str, suffix: str):
    path = Path(path_str)
    return sorted(
        p.name[:-len(suffix)]
        for p in path.glob(f"*{suffix}")
        if p.is_file() and p.stat().st_size > 0
    )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hmm_dir", required=True)
    p.add_argument("--tbl_dir", required=True)
    p.add_argument("--marker_dir", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    hmms = stems_from_dir(args.hmm_dir, ".hmm")
    tbls = stems_from_dir(args.tbl_dir, ".tbl")
    markers = stems_from_dir(args.marker_dir, ".fa")

    set_h = set(hmms)
    set_t = set(tbls)
    set_m = set(markers)
    common = sorted(set_h & set_t & set_m)

    with open(args.output, "w") as fh:
        fh.write("category\tog_id\n")
        for x in common:
            fh.write(f"common\t{x}\n")
        for x in sorted(set_h - set_t):
            fh.write(f"missing_tbl\t{x}\n")
        for x in sorted(set_h - set_m):
            fh.write(f"missing_marker\t{x}\n")
        for x in sorted(set_t - set_h):
            fh.write(f"missing_hmm\t{x}\n")
        for x in sorted(set_m - set_h):
            fh.write(f"missing_hmm_from_marker\t{x}\n")

    if not common:
        raise SystemExit("[ERROR] No OG IDs are shared across HMMs, tbls, and marker FASTAs.")

    if set_h != set_t or set_h != set_m:
        msg = (
            "[ERROR] Export input contract mismatch across HMMs / tbl / marker FASTAs.\n"
            f"        HMM count   : {len(set_h)}\n"
            f"        tbl count   : {len(set_t)}\n"
            f"        marker count: {len(set_m)}\n"
            f"        shared      : {len(common)}"
        )
        raise SystemExit(msg)

    print(f"  Export input audit passed: {len(common)} OG IDs shared across HMMs, tbls, and marker FASTAs.")

if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/audit_export_inputs.py"

cat > "${SCRIPTS_DIR}/filter_markers_with_outgroups.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
import math

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
    if cell is None:
        return 0
    value = str(cell).strip()
    if value == "" or value.lower() in {"nan", "none"}:
        return 0
    return len([x for x in value.split(",") if x.strip()])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--orthogroups_tsv", required=True)
    p.add_argument("--selected_ogs", required=True)
    p.add_argument("--retained_out", required=True)
    p.add_argument("--report_out", required=True)
    p.add_argument("--max_outgroup_presence_frac", type=float, default=0.40)
    p.add_argument("--max_outgroup_singlecopy_frac", type=float, default=0.50)
    args = p.parse_args()

    selected = set()
    with open(args.selected_ogs) as fh:
        for line in fh:
            og = line.strip()
            if og:
                selected.add(og)

    if not selected:
        raise SystemExit("[ERROR] Selected OG list is empty before outgroup filtering.")

    with open(args.orthogroups_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if not fieldnames:
            raise SystemExit("[ERROR] Orthogroups.tsv has no header.")
        og_field = fieldnames[0]

        cestoda_cols = [c for c in fieldnames if c in CESTODA]
        outgroup_cols = [c for c in fieldnames if c in OUTGROUPS]
        if not cestoda_cols:
            raise SystemExit("[ERROR] No Cestoda columns detected in Orthogroups.tsv")
        if not outgroup_cols:
            raise SystemExit("[ERROR] No outgroup columns detected in Orthogroups.tsv")

        n_out = len(outgroup_cols)
        max_out_presence_abs = math.floor(n_out * args.max_outgroup_presence_frac + 1e-9)
        rows = []
        kept = []

        for row in reader:
            og = row[og_field].strip()
            if og not in selected:
                continue

            out_present = 0
            out_singlecopy = 0
            out_total_copies = 0
            for sp in outgroup_cols:
                c = count_copies(row.get(sp))
                out_total_copies += c
                if c >= 1:
                    out_present += 1
                if c == 1:
                    out_singlecopy += 1

            out_presence_frac = out_present / n_out if n_out else 0.0
            out_singlecopy_frac = (out_singlecopy / out_present) if out_present else 0.0
            keep = (
                out_present <= max_out_presence_abs and
                out_singlecopy_frac <= args.max_outgroup_singlecopy_frac
            )

            rows.append({
                "orthogroup_id": og,
                "n_outgroups": n_out,
                "n_outgroup_present": out_present,
                "outgroup_presence_fraction": f"{out_presence_frac:.4f}",
                "n_outgroup_singlecopy": out_singlecopy,
                "outgroup_singlecopy_fraction_among_present": f"{out_singlecopy_frac:.4f}",
                "outgroup_total_copies": out_total_copies,
                "decision": "keep" if keep else "drop",
                "reason": (
                    "passes_outgroup_filter"
                    if keep else
                    f"outgroup_present>{max_out_presence_abs} or singlecopy_frac>{args.max_outgroup_singlecopy_frac:.2f}"
                ),
            })

            if keep:
                kept.append(og)

    if not kept:
        raise SystemExit("[ERROR] Outgroup filter removed all markers. Thresholds are too strict.")

    with open(args.retained_out, "w") as fh:
        for og in sorted(kept):
            fh.write(f"{og}\n")

    with open(args.report_out, "w", newline="") as fh:
        fieldnames = [
            "orthogroup_id",
            "n_outgroups",
            "n_outgroup_present",
            "outgroup_presence_fraction",
            "n_outgroup_singlecopy",
            "outgroup_singlecopy_fraction_among_present",
            "outgroup_total_copies",
            "decision",
            "reason",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"  Outgroup filter kept {len(kept)} of {len(selected)} candidate markers.")

if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/filter_markers_with_outgroups.py"

cat > "${SCRIPTS_DIR}/build_busco_id_map.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import csv

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--selected_ogs", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--prefix", type=int, default=1000000000)
    args = p.parse_args()

    with open(args.selected_ogs) as fh:
        og_ids = sorted({line.strip() for line in fh if line.strip()})
    if not og_ids:
        raise SystemExit("[ERROR] No retained OG IDs for BUSCO ID mapping.")

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["orthogroup_id", "busco_id"])
        for i, og in enumerate(og_ids, start=1):
            writer.writerow([og, str(args.prefix + i)])

    print(f"  Written BUSCO ID map: {args.output} ({len(og_ids)} entries)")

if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/build_busco_id_map.py"

cat > "${SCRIPTS_DIR}/compute_score_cutoffs_from_map.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def parse_tblout(path: Path):
    scores = []
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

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--id_map", required=True)
    p.add_argument("--tbl_dir", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--fraction", type=float, default=0.90)
    args = p.parse_args()

    mapping = []
    with open(args.id_map, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping.append((row["orthogroup_id"], row["busco_id"]))
    if not mapping:
        raise SystemExit("[ERROR] Empty BUSCO ID map.")

    cutoffs = {}
    for og_id, busco_id in mapping:
        tbl = Path(args.tbl_dir) / f"{og_id}.tbl"
        if not tbl.exists():
            raise SystemExit(f"[ERROR] Missing hmmsearch table for {og_id}: {tbl}")
        scores = parse_tblout(tbl)
        cutoffs[busco_id] = 999.0 if not scores else round(min(scores) * args.fraction, 2)

    with open(args.output, "w") as out:
        for busco_id in sorted(cutoffs, key=lambda x: int(x)):
            out.write(f"{busco_id}\t{cutoffs[busco_id]}\n")

    print(f"  Written score cutoffs: {len(cutoffs)}")

if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/compute_score_cutoffs_from_map.py"

cat > "${SCRIPTS_DIR}/rewrite_export_hmms.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

def rewrite_hmm_name(src: Path, dst: Path, new_name: str):
    seen_name = False
    with src.open() as infh, dst.open("w") as outfh:
        for line in infh:
            if line.startswith("NAME"):
                outfh.write(f"NAME  {new_name}\n")
                seen_name = True
            else:
                outfh.write(line)
        if not seen_name:
            raise SystemExit(f"[ERROR] HMM missing NAME field: {src}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--id_map", required=True)
    p.add_argument("--hmm_dir", required=True)
    p.add_argument("--out_dir", required=True)
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    with open(args.id_map, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            og = row["orthogroup_id"]
            busco = row["busco_id"]
            src = Path(args.hmm_dir) / f"{og}.hmm"
            dst = out_dir / f"{busco}.hmm"
            if not src.exists():
                raise SystemExit(f"[ERROR] Missing source HMM: {src}")
            rewrite_hmm_name(src, dst, busco)
            n += 1

    print(f"  Rewritten/exported HMMs: {n}")

if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/rewrite_export_hmms.py"

cat > "${SCRIPTS_DIR}/build_links_refseq_ancestral.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
import gzip
import subprocess
from pathlib import Path

def parse_fasta(path: Path):
    name = None
    seq_parts = []
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

def load_id_map(path: Path):
    ordered = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ordered.append((row["orthogroup_id"], row["busco_id"]))
    return ordered

def load_metadata(meta_path: Path):
    descriptions = {}
    if not meta_path.exists():
        return descriptions
    with meta_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if not fieldnames:
            return descriptions
        og_field = None
        for candidate in ("orthogroup_id", "og_id", "orthogroup", "busco_id"):
            if candidate in fieldnames:
                og_field = candidate
                break
        desc_field = None
        for candidate in ("conservation", "annotation", "description", "label"):
            if candidate in fieldnames:
                desc_field = candidate
                break
        if og_field is None:
            return descriptions
        for row in reader:
            og = row.get(og_field, "").strip()
            if og:
                descriptions[og] = row.get(desc_field, "").strip() if desc_field else "custom_cestoda_marker"
    return descriptions

def load_outgroup_decisions(path: Path):
    info = {}
    if not path.exists():
        return info
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            og = row.get("orthogroup_id", "").strip()
            if og:
                info[og] = row
    return info

def best_representative(marker_fa: Path):
    best_header = None
    best_seq = ""
    best_len = -1
    for header, seq in parse_fasta(marker_fa):
        if seq and len(seq) > best_len:
            best_header = header
            best_seq = seq
            best_len = len(seq)
    return best_header, best_seq

def emit_hmm_consensus(hmm_path: Path, busco_id: str):
    proc = subprocess.run(["hmmemit", "-c", str(hmm_path)], capture_output=True, text=True, check=True)
    out_lines = []
    for line in proc.stdout.splitlines():
        out_lines.append(f">{busco_id}" if line.startswith(">") else line)
    return "\n".join(out_lines) + "\n"

def emit_hmm_variants(hmm_path: Path, busco_id: str, n_variants: int = 5):
    proc = subprocess.run(["hmmemit", "-N", str(n_variants), str(hmm_path)], capture_output=True, text=True, check=True)
    out_lines = []
    idx = 0
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            idx += 1
            out_lines.append(f">{busco_id}_{idx}")
        else:
            out_lines.append(line)
    return "\n".join(out_lines) + "\n"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--id_map", required=True)
    p.add_argument("--marker_dir", required=True)
    p.add_argument("--hmm_dir", required=True)
    p.add_argument("--meta", required=True)
    p.add_argument("--outgroup_report", required=True)
    p.add_argument("--links_out", required=True)
    p.add_argument("--refseq_out", required=True)
    p.add_argument("--provenance_out", required=True)
    p.add_argument("--ancestral_out", required=True)
    p.add_argument("--ancestral_variants_out", required=True)
    args = p.parse_args()

    ordered = load_id_map(Path(args.id_map))
    descriptions = load_metadata(Path(args.meta))
    outgroup_info = load_outgroup_decisions(Path(args.outgroup_report))

    with open(args.links_out, "w") as links_fh,          gzip.open(args.refseq_out, "wt") as ref_fh,          open(args.provenance_out, "w") as prov_fh,          open(args.ancestral_out, "w") as anc_fh,          open(args.ancestral_variants_out, "w") as anc_var_fh:

        prov_fh.write(
            "busco_id\torthogroup_id\trepresentative_header\trepresentative_length\t"
            "outgroup_presence_fraction\toutgroup_singlecopy_fraction_among_present\n"
        )

        for og_id, busco_id in ordered:
            marker_fa = Path(args.marker_dir) / f"{og_id}.fa"
            hmm = Path(args.hmm_dir) / f"{busco_id}.hmm"
            if not marker_fa.exists():
                raise SystemExit(f"[ERROR] Missing marker FASTA for {og_id}: {marker_fa}")
            if not hmm.exists():
                raise SystemExit(f"[ERROR] Missing rewritten exported HMM for {busco_id}: {hmm}")

            desc = descriptions.get(og_id, "custom_cestoda_marker")
            out_meta = outgroup_info.get(og_id, {})
            out_pres = out_meta.get("outgroup_presence_fraction", "NA")
            out_sc = out_meta.get("outgroup_singlecopy_fraction_among_present", "NA")
            link = f"https://www.orthodb.org/?query={og_id}"
            links_fh.write(f"{busco_id}\t{desc}\t{link}\n")

            rep_header, rep_seq = best_representative(marker_fa)
            if not rep_seq:
                raise SystemExit(f"[ERROR] No valid protein sequences found in {marker_fa}")

            desc_text = (rep_header or "representative").replace("\t", " ").replace("\n", " ")
            ref_fh.write(f">{busco_id} representative={desc_text}\n")
            for i in range(0, len(rep_seq), 80):
                ref_fh.write(rep_seq[i:i+80] + "\n")

            prov_fh.write(f"{busco_id}\t{og_id}\t{desc_text}\t{len(rep_seq)}\t{out_pres}\t{out_sc}\n")
            anc_fh.write(emit_hmm_consensus(hmm, busco_id))
            anc_var_fh.write(emit_hmm_variants(hmm, busco_id, 5))

    print(f"  Written links      : {args.links_out}")
    print(f"  Written refseq DB  : {args.refseq_out}")
    print(f"  Written provenance : {args.provenance_out}")
    print(f"  Written ancestral  : {args.ancestral_out}")
    print(f"  Written variants   : {args.ancestral_variants_out}")

if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/build_links_refseq_ancestral.py"

cat > "${SCRIPTS_DIR}/final_lineage_contract_audit.py" << 'PYEOF'
#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path

def load_map(path: Path):
    ordered = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ordered.append((row["orthogroup_id"], row["busco_id"]))
    return ordered

def load_scores(path: Path):
    ids = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                ids.append(line.split("\t")[0])
    return ids

def load_ogs(path: Path):
    with path.open() as fh:
        return [line.strip() for line in fh if line.strip()]

def load_links(path: Path):
    ids = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                ids.append(line.split("\t")[0])
    return ids

def load_refseq_ids(path: Path):
    ids = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def hmm_name(path: Path):
    with path.open() as fh:
        for line in fh:
            if line.startswith("NAME"):
                return line.split(None, 1)[1].strip()
    return None

def count_headers(path: Path):
    c = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                c += 1
    return c

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--final_dir", required=True)
    p.add_argument("--id_map", required=True)
    args = p.parse_args()

    final_dir = Path(args.final_dir)
    mapping = load_map(Path(args.id_map))
    if not mapping:
        raise SystemExit("[ERROR] Empty final BUSCO ID map.")
    busco_ids = [b for _, b in mapping]
    expected = set(busco_ids)

    hmms_dir = final_dir / "hmms"
    if set(sorted(p.stem for p in hmms_dir.glob("*.hmm"))) != expected:
        raise SystemExit("[ERROR] Final HMM filenames do not match BUSCO ID map.")
    for busco_id in busco_ids:
        internal_name = hmm_name(hmms_dir / f"{busco_id}.hmm")
        if internal_name != busco_id:
            raise SystemExit(f"[ERROR] HMM NAME mismatch for {busco_id}: NAME={internal_name}")
    if set(load_scores(final_dir / "scores_cutoff")) != expected:
        raise SystemExit("[ERROR] scores_cutoff IDs do not match BUSCO ID map.")
    if set(load_ogs(final_dir / "info" / "ogs.id.info")) != expected:
        raise SystemExit("[ERROR] info/ogs.id.info IDs do not match BUSCO ID map.")
    if set(load_links(final_dir / "links_to_ODB12.txt")) != expected:
        raise SystemExit("[ERROR] links_to_ODB12.txt IDs do not match BUSCO ID map.")
    if set(load_refseq_ids(final_dir / "refseq_db.faa.gz")) != expected:
        raise SystemExit("[ERROR] refseq_db.faa.gz FASTA IDs do not match BUSCO ID map.")
    if count_headers(final_dir / "ancestral") != len(expected):
        raise SystemExit("[ERROR] ancestral sequence count mismatch.")
    if not (final_dir / "dataset.cfg").exists():
        raise SystemExit("[ERROR] Missing dataset.cfg")
    print(f"  Contract audit passed for {len(expected)} BUSCO markers.")

if __name__ == "__main__":
    main()
PYEOF
chmod +x "${SCRIPTS_DIR}/final_lineage_contract_audit.py"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Export input audit + BUSCO ID map"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

EXPORT_AUDIT_TSV="${FINAL_ODB_DIR}/info/export_input_audit.tsv"
python "${SCRIPTS_DIR}/audit_export_inputs.py"     --hmm_dir "${HMM_DIR}"     --tbl_dir "${HMMSEARCH_DIR}"     --marker_dir "${MARKER_DIR}"     --output "${EXPORT_AUDIT_TSV}"

SELECTED_OGS_ALL="${FINAL_ODB_DIR}/info/selected_ogs_pre_outgroup_filter.txt"
find "${MARKER_DIR}" -maxdepth 1 -type f -name "OG*.fa" -printf '%f\n' | sed 's/\.fa$//' | sort -u > "${SELECTED_OGS_ALL}"

OUTGROUP_FILTER_REPORT="${FINAL_ODB_DIR}/info/outgroup_filter_report.tsv"
SELECTED_OGS_FINAL="${FINAL_ODB_DIR}/info/selected_ogs_final.txt"

python "${SCRIPTS_DIR}/filter_markers_with_outgroups.py"     --orthogroups_tsv "${ORTHOGROUPS_TSV}"     --selected_ogs "${SELECTED_OGS_ALL}"     --retained_out "${SELECTED_OGS_FINAL}"     --report_out "${OUTGROUP_FILTER_REPORT}"     --max_outgroup_presence_frac 0.40     --max_outgroup_singlecopy_frac 0.50

N_SELECTED_PRE=$(wc -l < "${SELECTED_OGS_ALL}" | tr -d ' ')
N_SELECTED_POST=$(wc -l < "${SELECTED_OGS_FINAL}" | tr -d ' ')
N_FILTERED_OUT=$((N_SELECTED_PRE - N_SELECTED_POST))
echo "  Candidate markers before outgroup filter : ${N_SELECTED_PRE}"
echo "  Retained after outgroup filter           : ${N_SELECTED_POST}"
echo "  Removed by outgroup filter               : ${N_FILTERED_OUT}"

BUSCO_ID_MAP="${FINAL_ODB_DIR}/info/busco_id_map.tsv"
python "${SCRIPTS_DIR}/build_busco_id_map.py"     --selected_ogs "${SELECTED_OGS_FINAL}"     --output "${BUSCO_ID_MAP}"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 2 — Score cutoff computation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/compute_score_cutoffs_from_map.py"     --id_map "${BUSCO_ID_MAP}"     --tbl_dir "${HMMSEARCH_DIR}"     --output "${FINAL_ODB_DIR}/scores_cutoff"     --fraction "${SCORE_CUTOFF_FRACTION}"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 3 — Final HMM export"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/rewrite_export_hmms.py"     --id_map "${BUSCO_ID_MAP}"     --hmm_dir "${HMM_DIR}"     --out_dir "${FINAL_ODB_DIR}/hmms"

find "${FINAL_ODB_DIR}/prfl" -type f -delete || true
tail -n +2 "${BUSCO_ID_MAP}" | awk -F'\t' '{print $2}' | while read -r BUSCO_ID; do
    : > "${FINAL_ODB_DIR}/prfl/${BUSCO_ID}.prfl"
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 4 — Lineage sequence assets"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

META="${MARKER_DIR}/marker_metadata.tsv"
LINKS="${FINAL_ODB_DIR}/links_to_ODB12.txt"
REFSEQ_DB="${FINAL_ODB_DIR}/refseq_db.faa.gz"
PROVENANCE="${FINAL_ODB_DIR}/info/refseq_provenance.tsv"
ANCESTRAL="${FINAL_ODB_DIR}/ancestral"
ANCESTRAL_VARS="${FINAL_ODB_DIR}/ancestral_variants"

python "${SCRIPTS_DIR}/build_links_refseq_ancestral.py"     --id_map "${BUSCO_ID_MAP}"     --marker_dir "${MARKER_DIR}"     --hmm_dir "${FINAL_ODB_DIR}/hmms"     --meta "${META}"     --outgroup_report "${OUTGROUP_FILTER_REPORT}"     --links_out "${LINKS}"     --refseq_out "${REFSEQ_DB}"     --provenance_out "${PROVENANCE}"     --ancestral_out "${ANCESTRAL}"     --ancestral_variants_out "${ANCESTRAL_VARS}"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 5 — Dataset metadata"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

N_BUSCOS=$(tail -n +2 "${BUSCO_ID_MAP}" | wc -l | tr -d ' ')
TODAY=$(date +%Y-%m-%d)

cat > "${FINAL_ODB_DIR}/dataset.cfg" << EOF
name=${FINAL_DATASET_NAME}
domain=${DATASET_DOMAIN}
creation_date=${TODAY}
number_of_BUSCOs=${N_BUSCOS}
number_of_species=11
EOF

awk -F'\t' 'NR>1 && $2 != "" {print $2}' "${BUSCO_ID_MAP}" | sort -n > "${FINAL_ODB_DIR}/info/ogs.id.info"

cat > "${FINAL_ODB_DIR}/info/species.info" << 'EOF'
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
EOF

cp "${OUTGROUP_FILTER_REPORT}" "${FINAL_ODB_DIR}/info/outgroup_filter_report.tsv"
cp "${SELECTED_OGS_FINAL}" "${FINAL_ODB_DIR}/info/selected_ogs_final.txt"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 6 — Final lineage contract audit"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "${SCRIPTS_DIR}/final_lineage_contract_audit.py"     --final_dir "${FINAL_ODB_DIR}"     --id_map "${BUSCO_ID_MAP}"

echo ""
echo "════════════════════════════════════════════════════════"
echo "  Custom Cestoda BUSCO lineage export complete"
echo "════════════════════════════════════════════════════════"
echo ""
echo "  Dataset location : ${FINAL_ODB_DIR}/"
echo "  BUSCO markers    : ${N_BUSCOS}"
echo "  ID map           : ${BUSCO_ID_MAP}"
echo "  Outgroup report  : ${OUTGROUP_FILTER_REPORT}"
echo ""
echo "  Next step: bash 04_test_cestoda_lineage.sh"
