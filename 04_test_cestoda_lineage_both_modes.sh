#!/usr/bin/env bash
# =============================================================================
# 04_test_cestoda_lineage_both_modes.sh
# -----------------------------------------------------------------------------
# Validate the custom Cestoda BUSCO lineage in BOTH protein mode and genome mode
# for the SAME species used during lineage construction.
#
# Inputs:
#   - proteins: data/03_filtered_proteomes/<tag>.fa
#   - genomes : data/00_validation_genomes/<tag>.fa.gz
#
# Species list:
#   - logs/testing_species_manifest.tsv when present
#   - otherwise logs/species_manifest.tsv
#   - otherwise filtered proteome filenames
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

VALIDATION_GENOMES_DEFAULT="${PROJECT_ROOT}/data/00_validation_genomes"
GENOME_DIR="${GENOME_DIR:-${VALIDATION_GENOMES_DEFAULT}}"

mkdir -p "${VALIDATION_DIR}" "${LOG_DIR}"

LOG="${LOG_DIR}/04_test_cestoda_lineage_both_modes.log"
exec > >(tee -a "${LOG}") 2>&1

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
    BUSCO_CPUS=32
fi
if [[ "${BUSCO_CPUS}" -gt "${AVAILABLE_CORES}" ]]; then
    BUSCO_CPUS="${AVAILABLE_CORES}"
fi

DEFAULT_DATASET_PATH=""
if [[ -d "${PROJECT_ROOT}/${DATASET_NAME}_odb12" ]]; then
    DEFAULT_DATASET_PATH="${PROJECT_ROOT}/${DATASET_NAME}_odb12"
elif [[ -d "${ODB_DIR}" ]]; then
    DEFAULT_DATASET_PATH="${ODB_DIR}"
fi
DATASET_PATH="${DATASET_PATH:-${DEFAULT_DATASET_PATH}}"

if [[ -z "${DATASET_PATH}" ]] || [[ ! -d "${DATASET_PATH}" ]]; then
    echo "[ERROR] Validation dataset directory not found: ${DATASET_PATH}"
    exit 1
fi

DATASET_BASENAME="$(basename "${DATASET_PATH}")"
SUMMARY_TSV="${VALIDATION_DIR}/validation_both_modes_summary.tsv"

echo ""
echo "============================================================"
echo "  04_test_cestoda_lineage_both_modes.sh — $(date)"
echo "============================================================"
echo ""
echo "[INFO] Dataset         : ${DATASET_PATH}"
echo "[INFO] Validation root : ${VALIDATION_DIR}"
echo "[INFO] Genome dir      : ${GENOME_DIR}"
echo "[INFO] CPU cores       : ${AVAILABLE_CORES}"
echo "[INFO] BUSCO CPUs      : ${BUSCO_CPUS}"
echo ""

extract_summary_from_log() {
    local RUN_LOG="$1"
    awk '
        /\|C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+\s*\|/ {
            line=$0
            gsub(/^[[:space:]]*\|/, "", line)
            gsub(/\|[[:space:]]*$/, "", line)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
            print line
            exit
        }
    ' "${RUN_LOG}"
}

extract_summary_from_txt() {
    local SUMMARY_TXT="$1"
    awk '
        /\#.*C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+/ {
            line=$0
            sub(/^#*[[:space:]]*/, "", line)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
            print line
            exit
        }
    ' "${SUMMARY_TXT}"
}

extract_summary_from_json() {
    local SUMMARY_JSON="$1"
    python - "$SUMMARY_JSON" << 'PYEOF'
import json
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
pat = re.compile(r'C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+')

def walk(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from walk(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk(v)
    elif isinstance(obj, str):
        yield obj

with path.open() as fh:
    data = json.load(fh)

for value in walk(data):
    m = pat.search(value)
    if m:
        print(m.group(0))
        sys.exit(0)

sys.exit(1)
PYEOF
}

get_busco_summary() {
    local RUN_DIR="$1"
    local RUN_LOG="$2"
    local SUMMARY=""

    SUMMARY="$(extract_summary_from_log "${RUN_LOG}" || true)"
    if [[ -n "${SUMMARY}" ]]; then
        printf '%s\n' "${SUMMARY}"
        return 0
    fi

    local TXT
    TXT="$(find "${RUN_DIR}" -maxdepth 2 -type f -name 'short_summary*.txt' | sort | head -n 1 || true)"
    if [[ -n "${TXT}" ]] && [[ -f "${TXT}" ]]; then
        SUMMARY="$(extract_summary_from_txt "${TXT}" || true)"
        if [[ -n "${SUMMARY}" ]]; then
            printf '%s\n' "${SUMMARY}"
            return 0
        fi
    fi

    local JSON
    JSON="$(find "${RUN_DIR}" -maxdepth 2 -type f -name 'short_summary*.json' | sort | head -n 1 || true)"
    if [[ -n "${JSON}" ]] && [[ -f "${JSON}" ]]; then
        SUMMARY="$(extract_summary_from_json "${JSON}" || true)"
        if [[ -n "${SUMMARY}" ]]; then
            printf '%s\n' "${SUMMARY}"
            return 0
        fi
    fi

    return 1
}

find_genome_for_species() {
    local SPECIES="$1"
    local -a EXTS=("fa" "fna" "fasta" "fa.gz" "fna.gz" "fasta.gz")
    local EXT MATCH

    [[ -d "${GENOME_DIR}" ]] || return 1

    for EXT in "${EXTS[@]}"; do
        while IFS= read -r MATCH; do
            [[ -n "${MATCH}" ]] || continue
            if [[ -s "${MATCH}" ]]; then
                printf '%s\n' "${MATCH}"
                return 0
            fi
        done < <(find "${GENOME_DIR}" -maxdepth 3 -type f \
            \( -iname "${SPECIES}.${EXT}" -o -iname "${SPECIES}_*.${EXT}" -o -iname "*${SPECIES}*.${EXT}" \) \
            | sort)
    done

    return 1
}

run_busco_validation() {
    local SPECIES="$1"
    local MODE="$2"
    local INPUT="$3"

    local RUN_NAME="validate_${SPECIES}_${MODE}_${DATASET_BASENAME}"
    local OUT_DIR="${VALIDATION_DIR}/${RUN_NAME}"
    local RUN_LOG="${VALIDATION_DIR}/${RUN_NAME}.console.log"

    rm -rf "${OUT_DIR}"
    rm -f "${RUN_LOG}"

    echo "  [RUN] ${RUN_NAME} (${MODE})"
    if ! busco \
        --cpu "${BUSCO_CPUS}" \
        -i "${INPUT}" \
        -m "${MODE}" \
        -l "${DATASET_PATH}" \
        -o "${RUN_NAME}" \
        --out_path "${VALIDATION_DIR}" \
        --offline -f 2>&1 | tee "${RUN_LOG}"; then
        echo "[ERROR] BUSCO failed for ${SPECIES} (${MODE})"
        return 1
    fi

    if grep -Fq "No jobs to run on hmmsearch" "${RUN_LOG}"; then
        echo "[ERROR] Historical failure detected for ${SPECIES} (${MODE}): No jobs to run on hmmsearch"
        return 1
    fi

    local RUN_SUBDIR
    RUN_SUBDIR="$(find "${OUT_DIR}" -maxdepth 1 -type d -name 'run_*' | sort | head -n 1 || true)"
    if [[ -z "${RUN_SUBDIR}" ]] || [[ ! -d "${RUN_SUBDIR}" ]]; then
        echo "[ERROR] BUSCO run directory not found under ${OUT_DIR}"
        return 1
    fi

    local SUMMARY
    if ! SUMMARY="$(get_busco_summary "${RUN_SUBDIR}" "${RUN_LOG}")"; then
        echo "[ERROR] Could not parse BUSCO summary from ${RUN_SUBDIR}"
        return 1
    fi

    echo "  [OK] ${SPECIES} (${MODE}) -> ${SUMMARY}"
    printf "%s\t%s\t%s\t%s\t%s\n" \
        "${SPECIES}" "${MODE}" "${INPUT}" "${RUN_NAME}" "${SUMMARY}" \
        >> "${SUMMARY_TSV}"
}

build_species_list() {
    local TMP_LIST="$1"
    local TEST_MANIFEST="${LOG_DIR}/testing_species_manifest.tsv"
    local MAIN_MANIFEST="${LOG_DIR}/species_manifest.tsv"

    : > "${TMP_LIST}"

    if [[ -s "${TEST_MANIFEST}" ]]; then
        awk -F'\t' 'NR>1 && $1 != "" {print $1}' "${TEST_MANIFEST}" | sort -u > "${TMP_LIST}"
    elif [[ -s "${MAIN_MANIFEST}" ]]; then
        awk -F'\t' 'NR>1 && $1 != "" {print $1}' "${MAIN_MANIFEST}" | sort -u > "${TMP_LIST}"
    else
        find "${FILTERED_PROTEOMES}" -maxdepth 1 -type f -name '*.fa' -printf '%f\n' \
            | sed 's/\.fa$//' | sort -u > "${TMP_LIST}"
    fi
}

printf "species\tmode\tinput_path\trun_name\tsummary\n" > "${SUMMARY_TSV}"
TMP_SPECIES="$(mktemp)"
trap 'rm -f "${TMP_SPECIES}"' EXIT
build_species_list "${TMP_SPECIES}"

TOTAL_SPECIES="$(wc -l < "${TMP_SPECIES}" | tr -d ' ')"
MISSING_PROTEOMES=0
MISSING_GENOMES=0

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Protein- and genome-mode validation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "[INFO] Species to validate: ${TOTAL_SPECIES}"
echo ""

while IFS= read -r SPECIES; do
    [[ -n "${SPECIES}" ]] || continue

    PROTEOME="${FILTERED_PROTEOMES}/${SPECIES}.fa"
    GENOME="$(find_genome_for_species "${SPECIES}" || true)"

    if [[ ! -s "${PROTEOME}" ]]; then
        echo "  [WARN] Missing proteome for ${SPECIES}: ${PROTEOME}"
        printf "%s\tprotein\t%s\t%s\t%s\n" \
            "${SPECIES}" "${PROTEOME}" "NA" "SKIPPED_missing_proteome" \
            >> "${SUMMARY_TSV}"
        MISSING_PROTEOMES=$((MISSING_PROTEOMES + 1))
    else
        run_busco_validation "${SPECIES}" "protein" "${PROTEOME}"
    fi

    if [[ -z "${GENOME}" ]] || [[ ! -s "${GENOME}" ]]; then
        echo "  [WARN] Missing genome for ${SPECIES}"
        printf "%s\tgenome\t%s\t%s\t%s\n" \
            "${SPECIES}" "${GENOME:-NA}" "NA" "SKIPPED_missing_genome" \
            >> "${SUMMARY_TSV}"
        MISSING_GENOMES=$((MISSING_GENOMES + 1))
    else
        run_busco_validation "${SPECIES}" "genome" "${GENOME}"
    fi

    echo ""
done < "${TMP_SPECIES}"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Validation complete"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Summary TSV        : ${SUMMARY_TSV}"
echo "  Missing proteomes  : ${MISSING_PROTEOMES}"
echo "  Missing genomes    : ${MISSING_GENOMES}"
echo ""

python - "${SUMMARY_TSV}" << 'PYEOF'
import sys
from pathlib import Path

tsv = Path(sys.argv[1])
print("  Parsed validation summaries:")
print("  ----------------------------")
with tsv.open() as fh:
    next(fh, None)
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 5:
            continue
        species, mode, input_path, run_name, summary = parts
        print(f"  {species:<28} {mode:<7} {summary}")
PYEOF

echo ""
echo "[OK] Validation script finished."
