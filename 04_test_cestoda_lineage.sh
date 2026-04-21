#!/usr/bin/env bash
# =============================================================================
# 04_test_cestoda_lineage.sh
# -----------------------------------------------------------------------------
# Functional validation for the exported custom Cestoda BUSCO lineage.
#
# Goals:
#   1. Verify BUSCO can run end-to-end with the local lineage.
#   2. Fail loudly if BUSCO reports "No jobs to run on hmmsearch".
#   3. Record robust summaries without depending on one fragile BUSCO JSON/TXT
#      schema.
#
# Validation plan:
#   - Protein mode on the species listed in VALIDATION_SPECIES
#   - Optional genome-mode validation if a known reference genome file exists
#
# Usage:
#   bash 04_test_cestoda_lineage.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

if command -v nproc >/dev/null 2>&1; then
    AVAILABLE_CORES="$(nproc)"
else
    AVAILABLE_CORES="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi

if [[ -z "${AVAILABLE_CORES}" ]] || ! [[ "${AVAILABLE_CORES}" =~ ^[0-9]+$ ]] || [[ "${AVAILABLE_CORES}" -lt 1 ]]; then
    AVAILABLE_CORES=1
fi

BUSCO_CPUS="${BUSCO_CPUS:-32}"
if ! [[ "${BUSCO_CPUS}" =~ ^[0-9]+$ ]] || [[ "${BUSCO_CPUS}" -lt 1 ]]; then
    BUSCO_CPUS=1
fi
if [[ "${BUSCO_CPUS}" -gt "${AVAILABLE_CORES}" ]]; then
    BUSCO_CPUS="${AVAILABLE_CORES}"
fi

if [[ "${DATASET_NAME}" == *_odb12 ]]; then
    FINAL_DATASET_NAME="${DATASET_NAME}"
else
    FINAL_DATASET_NAME="${DATASET_NAME}_odb12"
fi
DATASET_DIR="${PROJECT_ROOT}/${FINAL_DATASET_NAME}"

mkdir -p "${VALIDATION_DIR}" "${LOG_DIR}"
LOG="${LOG_DIR}/04_validation.log"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  04_test_cestoda_lineage.sh — $(date)"
echo "============================================================"
echo ""
echo "[INFO] Dataset         : ${DATASET_DIR}"
echo "[INFO] Validation root : ${VALIDATION_DIR}"
echo "[INFO] CPU cores       : ${AVAILABLE_CORES}"
echo "[INFO] BUSCO CPUs      : ${BUSCO_CPUS}"

if [[ ! -d "${DATASET_DIR}" ]]; then
    echo "[ERROR] Dataset directory not found: ${DATASET_DIR}"
    echo "        Run 03_build_odb.sh first."
    exit 1
fi

for req in dataset.cfg hmms scores_cutoff links_to_ODB12.txt refseq_db.faa.gz ancestral ancestral_variants; do
    if [[ ! -e "${DATASET_DIR}/${req}" ]]; then
        echo "[ERROR] Missing required dataset asset: ${DATASET_DIR}/${req}"
        exit 1
    fi
done

SUMMARY_TSV="${VALIDATION_DIR}/validation_summary.tsv"
echo -e "run_name\tmode\tinput\tstatus_line\tcomplete_pct\tsingle_pct\tdup_pct\tfrag_pct\tmissing_pct\tn_buscos" > "${SUMMARY_TSV}"

extract_summary_line() {
    local run_dir="$1"
    local console_log="$2"
    local summary_line=""

    if [[ -s "${console_log}" ]]; then
        summary_line="$(grep -E '\|C:[0-9]+(\.[0-9]+)?%\[S:[0-9]+(\.[0-9]+)?%,D:[0-9]+(\.[0-9]+)?%\],F:[0-9]+(\.[0-9]+)?%,M:[0-9]+(\.[0-9]+)?%,n:[0-9]+' "${console_log}" | tail -n 1 | sed 's/^ *|//; s/ *|$//')"
    fi

    if [[ -z "${summary_line}" ]]; then
        local txt_file=""
        txt_file="$(find "${run_dir}" -maxdepth 2 -type f -name 'short_summary*.txt' | sort | head -n 1 || true)"
        if [[ -n "${txt_file}" && -s "${txt_file}" ]]; then
            summary_line="$(grep -E '^C:[0-9]+(\.[0-9]+)?%\[S:[0-9]+(\.[0-9]+)?%,D:[0-9]+(\.[0-9]+)?%\],F:[0-9]+(\.[0-9]+)?%,M:[0-9]+(\.[0-9]+)?%,n:[0-9]+' "${txt_file}" | head -n 1 || true)"
        fi
    fi

    if [[ -z "${summary_line}" ]]; then
        local json_file=""
        json_file="$(find "${run_dir}" -maxdepth 2 -type f -name 'short_summary*.json' | sort | head -n 1 || true)"
        if [[ -n "${json_file}" && -s "${json_file}" ]]; then
            summary_line="$(JSON_FILE="${json_file}" python - <<'PYEOF'
import json
import os
import re
from collections.abc import Mapping, Sequence

path = os.environ['JSON_FILE']
with open(path) as fh:
    obj = json.load(fh)

pattern = re.compile(r'^C:[0-9]+(?:\.[0-9]+)?%\[S:[0-9]+(?:\.[0-9]+)?%,D:[0-9]+(?:\.[0-9]+)?%\],F:[0-9]+(?:\.[0-9]+)?%,M:[0-9]+(?:\.[0-9]+)?%,n:[0-9]+$')

seen = set()
stack = [obj]
while stack:
    cur = stack.pop()
    ident = id(cur)
    if ident in seen:
        continue
    seen.add(ident)
    if isinstance(cur, str) and pattern.match(cur.strip()):
        print(cur.strip())
        raise SystemExit(0)
    if isinstance(cur, Mapping):
        for v in cur.values():
            stack.append(v)
    elif isinstance(cur, Sequence) and not isinstance(cur, (str, bytes, bytearray)):
        for v in cur:
            stack.append(v)

# fallback: reconstruct from common BUSCO json keys if present anywhere
fields = {
    'complete': None,
    'single': None,
    'duplicated': None,
    'fragmented': None,
    'missing': None,
    'n_markers': None,
}

stack = [obj]
seen.clear()
while stack:
    cur = stack.pop()
    ident = id(cur)
    if ident in seen:
        continue
    seen.add(ident)
    if isinstance(cur, Mapping):
        low = {str(k).lower(): v for k, v in cur.items()}
        for k in list(fields):
            if fields[k] is None and k in low:
                fields[k] = low[k]
        for alt_src, alt_dst in [
            ('complete_buscos','complete'),
            ('single_copy_buscos','single'),
            ('multi_copy_buscos','duplicated'),
            ('fragmented_buscos','fragmented'),
            ('missing_buscos','missing'),
            ('n','n_markers'),
            ('total_buscos','n_markers'),
            ('total_busco_groups','n_markers'),
        ]:
            if fields[alt_dst] is None and alt_src in low:
                fields[alt_dst] = low[alt_src]
        for v in cur.values():
            stack.append(v)
    elif isinstance(cur, Sequence) and not isinstance(cur, (str, bytes, bytearray)):
        for v in cur:
            stack.append(v)

if all(fields[k] is not None for k in fields):
    def fmt(x):
        try:
            return f"{float(x):.1f}"
        except Exception:
            return str(x)
    print(
        f"C:{fmt(fields['complete'])}%[S:{fmt(fields['single'])}%,D:{fmt(fields['duplicated'])}%],"
        f"F:{fmt(fields['fragmented'])}%,M:{fmt(fields['missing'])}%,n:{int(float(fields['n_markers']))}",
        sep=''
    )
PYEOF
)"
        fi
    fi

    printf '%s\n' "${summary_line}"
}

parse_summary_fields() {
    local summary_line="$1"
    SUMMARY_LINE="${summary_line}" python - <<'PYEOF'
import os, re, sys
s = os.environ['SUMMARY_LINE'].strip()
m = re.match(r'^C:([0-9]+(?:\.[0-9]+)?)%\[S:([0-9]+(?:\.[0-9]+)?)%,D:([0-9]+(?:\.[0-9]+)?)%\],F:([0-9]+(?:\.[0-9]+)?)%,M:([0-9]+(?:\.[0-9]+)?)%,n:([0-9]+)$', s)
if not m:
    sys.exit(1)
print("\t".join(m.groups()))
PYEOF
}

run_busco_validation() {
    local mode="$1"
    local input_path="$2"
    local run_name="$3"

    if [[ ! -s "${input_path}" ]]; then
        echo "  [WARN] Missing input for ${run_name}: ${input_path}"
        return 0
    fi

    local out_dir="${VALIDATION_DIR}/${run_name}"
    local console_log="${out_dir}.busco.console.log"
    rm -rf "${out_dir}" "${console_log}"

    echo "  [RUN] ${run_name} (${mode})"
    busco \
        --cpu "${BUSCO_CPUS}" \
        -i "${input_path}" \
        -m "${mode}" \
        -l "${DATASET_DIR}" \
        -o "${run_name}" \
        --out_path "${VALIDATION_DIR}" \
        --offline -f 2>&1 | tee "${console_log}"

    if grep -Fq 'No jobs to run on hmmsearch' "${console_log}"; then
        echo "[ERROR] BUSCO reported 'No jobs to run on hmmsearch' for ${run_name}" 
        exit 1
    fi

    local run_subdir=""
    run_subdir="$(find "${out_dir}" -maxdepth 1 -mindepth 1 -type d -name 'run_*' | sort | head -n 1 || true)"
    if [[ -z "${run_subdir}" ]]; then
        echo "[ERROR] BUSCO run directory not found under ${out_dir}"
        exit 1
    fi

    local summary_line=""
    summary_line="$(extract_summary_line "${run_subdir}" "${console_log}")"
    if [[ -z "${summary_line}" ]]; then
        echo "[ERROR] Could not parse BUSCO summary from ${run_subdir}"
        echo "        Console log inspected: ${console_log}"
        exit 1
    fi

    local parsed=""
    if ! parsed="$(parse_summary_fields "${summary_line}")"; then
        echo "[ERROR] Parsed summary line has unexpected format: ${summary_line}"
        exit 1
    fi

    local c s d f m n
    IFS=$'\t' read -r c s d f m n <<< "${parsed}"
    printf "  [OK]  %s\n" "${summary_line}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${run_name}" "${mode}" "${input_path}" "${summary_line}" "${c}" "${s}" "${d}" "${f}" "${m}" "${n}" \
        >> "${SUMMARY_TSV}"
}

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Proteins-mode validation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

for species in "${VALIDATION_SPECIES[@]}"; do
    run_busco_validation "protein" "${FILTERED_PROTEOMES}/${species}.fa" "validate_${species}_${FINAL_DATASET_NAME}"
done

# Optional genome-mode smoke test on a known local genome if present.
GENOME_INPUT="${PROJECT_ROOT}/echinococcus_canadensis.PRJEB8992.WBPS19.genomic_masked.fa.gz"
if [[ -s "${GENOME_INPUT}" ]]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Stage 2 — Genome-mode validation"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    run_busco_validation "genome" "${GENOME_INPUT}" "validate_echinococcus_canadensis_genome_${FINAL_DATASET_NAME}"
else
    echo ""
    echo "[INFO] Genome validation input not found, skipping: ${GENOME_INPUT}"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Validation complete — $(date)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Summary table: ${SUMMARY_TSV}"
cat "${SUMMARY_TSV}"
