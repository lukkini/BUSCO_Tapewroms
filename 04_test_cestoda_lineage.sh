#!/usr/bin/env bash
# =============================================================================
# 04_test_cestoda_lineage.sh
# -----------------------------------------------------------------------------
# Validate the exported custom Cestoda BUSCO lineage in protein mode across the
# configured training/reference species.
#
# Goals:
#   - run BUSCO fresh for each validation species
#   - fail loudly if BUSCO reports "No jobs to run on hmmsearch"
#   - parse BUSCO v6 summaries robustly from short_summary*.txt or JSON
#   - write a machine-readable validation summary TSV
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

mkdir -p "${VALIDATION_DIR}" "${LOG_DIR}" "${SCRIPTS_DIR}"

LOG="${LOG_DIR}/04_validation.log"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  04_test_cestoda_lineage.sh — $(date)"
echo "============================================================"
echo ""
echo "[INFO] Dataset         : ${FINAL_ODB_DIR}"
echo "[INFO] Validation root : ${VALIDATION_DIR}"
if command -v nproc >/dev/null 2>&1; then
    echo "[INFO] CPU cores       : $(nproc)"
fi
echo ""

if [[ ! -d "${FINAL_ODB_DIR}" ]]; then
    echo "[ERROR] Final lineage directory not found:"
    echo "        ${FINAL_ODB_DIR}"
    echo "        Run 03_build_odb.sh first."
    exit 1
fi

if [[ ! -s "${FINAL_ODB_DIR}/dataset.cfg" ]]; then
    echo "[ERROR] Missing dataset.cfg in ${FINAL_ODB_DIR}"
    exit 1
fi

if [[ ! -d "${FINAL_ODB_DIR}/hmms" ]]; then
    echo "[ERROR] Missing hmms/ in ${FINAL_ODB_DIR}"
    exit 1
fi

if [[ ! -s "${FINAL_ODB_DIR}/scores_cutoff" ]]; then
    echo "[ERROR] Missing scores_cutoff in ${FINAL_ODB_DIR}"
    exit 1
fi

find_run_dir() {
    local out_name="$1"
    local out_root="${VALIDATION_DIR}/${out_name}"
    local candidate

    if [[ -d "${out_root}/run_${FINAL_DATASET_NAME}" ]]; then
        echo "${out_root}/run_${FINAL_DATASET_NAME}"
        return 0
    fi

    candidate="$(find "${out_root}" -maxdepth 1 -mindepth 1 -type d -name 'run_*' | head -n 1 || true)"
    if [[ -n "${candidate}" ]]; then
        echo "${candidate}"
        return 0
    fi

    return 1
}

find_short_summary() {
    local run_dir="$1"
    local out_root="$2"
    local candidate=""

    candidate="$(find "${run_dir}" -maxdepth 1 -type f -name 'short_summary*.json' | sort | head -n 1 || true)"
    if [[ -n "${candidate}" ]]; then
        echo "${candidate}"
        return 0
    fi

    candidate="$(find "${run_dir}" -maxdepth 1 -type f -name 'short_summary*.txt' | sort | head -n 1 || true)"
    if [[ -n "${candidate}" ]]; then
        echo "${candidate}"
        return 0
    fi

    candidate="$(find "${out_root}" -maxdepth 2 -type f -name 'short_summary*.json' | sort | head -n 1 || true)"
    if [[ -n "${candidate}" ]]; then
        echo "${candidate}"
        return 0
    fi

    candidate="$(find "${out_root}" -maxdepth 2 -type f -name 'short_summary*.txt' | sort | head -n 1 || true)"
    if [[ -n "${candidate}" ]]; then
        echo "${candidate}"
        return 0
    fi

    return 1
}

parse_busco_summary() {
    local summary_path="$1"

    python - "$summary_path" << 'PYEOF'
import json
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])

def emit(vals):
    print("\t".join(str(v) for v in vals))

if path.suffix == ".json":
    data = json.loads(path.read_text())
    results = data.get("results", {})
    percentages = results.get("Complete percentage")
    single_pct = results.get("Single copy percentage")
    dup_pct = results.get("Multi copy percentage")
    frag_pct = results.get("Fragmented percentage")
    miss_pct = results.get("Missing percentage")
    n = results.get("n_markers")
    complete = results.get("Complete", {}).get("count")
    single = results.get("Single copy", {}).get("count")
    dup = results.get("Multi copy", {}).get("count")
    frag = results.get("Fragmented", {}).get("count")
    miss = results.get("Missing", {}).get("count")
    if None in (percentages, single_pct, dup_pct, frag_pct, miss_pct, n, complete, single, dup, frag, miss):
        raise SystemExit(2)
    emit([percentages, single_pct, dup_pct, frag_pct, miss_pct, n, complete, single, dup, frag, miss])
    raise SystemExit(0)

text = path.read_text()

line = None
for raw in text.splitlines():
    s = raw.strip()
    if s.startswith("C:") and ",n:" in s:
        line = s
        break

if line is None:
    raise SystemExit(2)

m = re.search(
    r"C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%\],"
    r"F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%,n:(?P<n>[0-9]+)",
    line,
)
if not m:
    raise SystemExit(2)

counts = {
    "Complete BUSCOs": None,
    "Complete and single-copy BUSCOs": None,
    "Complete and duplicated BUSCOs": None,
    "Fragmented BUSCOs": None,
    "Missing BUSCOs": None,
}

for raw in text.splitlines():
    s = raw.strip()
    for label in list(counts):
        if label in s:
            m2 = re.match(r"(?P<count>[0-9]+)\s+", s)
            if m2:
                counts[label] = int(m2.group("count"))

if any(v is None for v in counts.values()):
    raise SystemExit(2)

emit([
    m.group("C"),
    m.group("S"),
    m.group("D"),
    m.group("F"),
    m.group("M"),
    m.group("n"),
    counts["Complete BUSCOs"],
    counts["Complete and single-copy BUSCOs"],
    counts["Complete and duplicated BUSCOs"],
    counts["Fragmented BUSCOs"],
    counts["Missing BUSCOs"],
])
PYEOF
}

assert_no_hmmsearch_contract_failure() {
    local out_root="$1"
    if grep -R -F -q "No jobs to run on hmmsearch" "${out_root}" 2>/dev/null; then
        echo "[ERROR] BUSCO reported 'No jobs to run on hmmsearch' in:"
        echo "        ${out_root}"
        exit 1
    fi
}

SUMMARY_TSV="${VALIDATION_DIR}/validation_summary.tsv"
printf "species\tmode\tcomplete_pct\tsingle_pct\tduplicated_pct\tfragmented_pct\tmissing_pct\tn_markers\tcomplete\tsingle\tduplicated\tfragmented\tmissing\tsummary_file\n" > "${SUMMARY_TSV}"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Proteins-mode validation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

N_DONE=0

for SPECIES in "${VALIDATION_SPECIES[@]}"; do
    INPUT_FA="${FILTERED_PROTEOMES}/${SPECIES}.fa"
    OUT_NAME="validate_${SPECIES}_${FINAL_DATASET_NAME}"
    OUT_ROOT="${VALIDATION_DIR}/${OUT_NAME}"

    if [[ ! -s "${INPUT_FA}" ]]; then
        echo "[ERROR] Validation proteome missing: ${INPUT_FA}"
        exit 1
    fi

    echo "  [RUN] ${OUT_NAME} (proteins)"
    rm -rf "${OUT_ROOT}"

    busco         --cpu "${OF_THREADS}"         -i "${INPUT_FA}"         -m protein         -l "${FINAL_ODB_DIR}"         -o "${OUT_NAME}"         --out_path "${VALIDATION_DIR}"         --offline         -f

    assert_no_hmmsearch_contract_failure "${OUT_ROOT}"

    RUN_DIR="$(find_run_dir "${OUT_NAME}")" || {
        echo "[ERROR] Could not locate BUSCO run directory under ${OUT_ROOT}"
        exit 1
    }

    SUMMARY_FILE="$(find_short_summary "${RUN_DIR}" "${OUT_ROOT}")" || {
        echo "[ERROR] Could not locate BUSCO short summary under ${OUT_ROOT}"
        exit 1
    }

    if ! PARSED="$(parse_busco_summary "${SUMMARY_FILE}")"; then
        echo "[ERROR] Could not parse BUSCO summary from ${SUMMARY_FILE}"
        exit 1
    fi

    IFS=$'\t' read -r C_PCT S_PCT D_PCT F_PCT M_PCT N_MARKERS N_COMPLETE N_SINGLE N_DUP N_FRAG N_MISS <<< "${PARSED}"

    printf "%s\tprotein\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"         "${SPECIES}"         "${C_PCT}" "${S_PCT}" "${D_PCT}" "${F_PCT}" "${M_PCT}"         "${N_MARKERS}" "${N_COMPLETE}" "${N_SINGLE}" "${N_DUP}" "${N_FRAG}" "${N_MISS}"         "${SUMMARY_FILE}"         >> "${SUMMARY_TSV}"

    printf "        C:%s%%  S:%s%%  D:%s%%  F:%s%%  M:%s%%  n:%s\n"         "${C_PCT}" "${S_PCT}" "${D_PCT}" "${F_PCT}" "${M_PCT}" "${N_MARKERS}"

    N_DONE=$((N_DONE + 1))
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Validation complete — $(date)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Validation runs completed : ${N_DONE}"
echo "  Summary TSV               : ${SUMMARY_TSV}"
echo "  Log                       : ${LOG}"
echo ""
