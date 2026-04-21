#!/usr/bin/env bash
# =============================================================================
# 00_setup_environment.sh
# -----------------------------------------------------------------------------
# Creates the conda environment and all configuration files for the Cestoda
# BUSCO lineage pipeline.
#
# Assumptions:
#   - Anaconda is installed and initialised (conda command available)
#   - conda is configured with the libmamba solver (default in Anaconda ≥ 2022)
#
# Usage:
#   bash 00_setup_environment.sh
#
# After completion:
#   conda activate cestoda_busco_builder
#   bash 01_retrieve_data.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# STEP 1 — Write config/project.cfg
# =============================================================================

mkdir -p "${SCRIPT_DIR}/config"
mkdir -p "${SCRIPT_DIR}/logs"

cat > "${SCRIPT_DIR}/config/project.cfg" << 'CFGEOF'
# =============================================================================
# config/project.cfg — shared configuration for the Cestoda BUSCO pipeline.
# Sourced by every pipeline script.
# =============================================================================

CONDA_ENV="cestoda_busco_builder"

# Project root is the parent of this config/ directory
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---- Directory layout -------------------------------------------------------
RAW_PROTEOMES="${PROJECT_ROOT}/data/01_raw_proteomes"
CLEAN_PROTEOMES="${PROJECT_ROOT}/data/02_clean_proteomes"
FILTERED_PROTEOMES="${PROJECT_ROOT}/data/03_filtered_proteomes"
OG_DIR="${PROJECT_ROOT}/data/04_orthofinder"
MARKER_DIR="${PROJECT_ROOT}/data/05_cestoda_markers"
ALIGN_DIR="${PROJECT_ROOT}/data/06_alignments"
HMM_DIR="${PROJECT_ROOT}/data/07_hmm_profiles"
HMMSEARCH_DIR="${PROJECT_ROOT}/data/08_hmmsearch"
ODB_DIR="${PROJECT_ROOT}/cestoda_odb_custom"
VALIDATION_DIR="${PROJECT_ROOT}/validation"
LOG_DIR="${PROJECT_ROOT}/logs"
SCRIPTS_DIR="${PROJECT_ROOT}/scripts"

# ---- WormBase ParaSite release ----------------------------------------------
# Verify the current release at: https://parasite.wormbase.org/ftp.html
WBPS_RELEASE="WBPS19"
WBPS_BASE="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${WBPS_RELEASE}/species"

# ---- Compute resources ------------------------------------------------------
# Edit OF_THREADS to match your server. Check available CPUs with: nproc
OF_THREADS=32
OF_ANALYSIS_THREADS=8
PARALLEL_JOBS=8
THREADS_PER_JOB=4

# ---- Sequence quality-control thresholds ------------------------------------
MIN_PROTEIN_LENGTH=30
MAX_X_FRACTION=0.10

# ---- Marker-selection thresholds --------------------------------------------
MIN_CESTODA_PRESENCE=0.80
MIN_SINGLECOPY_FRAC=0.90
OUTGROUP_CONSERVATION_MIN=3

# ---- Score cutoff -----------------------------------------------------------
SCORE_CUTOFF_FRACTION=0.90

# ---- Download settings ------------------------------------------------------
WGET_RETRIES=5
WGET_TIMEOUT=60

# ---- Dataset metadata -------------------------------------------------------
DATASET_NAME="cestoda_odb_custom"
DATASET_LINEAGE="Cestoda"
DATASET_DOMAIN="eukaryota"

# ---- Validation species -----------------------------------------------------
VALIDATION_SPECIES=(
    "echinococcus_multilocularis"
    "hymenolepis_microstoma"
    "taenia_multiceps"
)

# ---- Cestoda species tags (for refseq_db construction) ----------------------
CESTODA_TAGS=(
    echinococcus_canadensis_g7
    echinococcus_granulosus_g1
    echinococcus_multilocularis
    hymenolepis_diminuta
    hymenolepis_microstoma
    hymenolepis_nana
    mesocestoides_corti
    taenia_asiatica
    taenia_multiceps
    taenia_saginata
    taenia_solium
)
CFGEOF

echo "[OK] config/project.cfg written."

# =============================================================================
# STEP 2 — Write config/environment.yml
# =============================================================================
# Version pinning notes:
#   - busco=6.0.0 + sepp=4.5.5  : sepp 4.5.6 breaks BUSCO v6; keep both pinned
#   - orthofinder                : NOT pinned — bioconda's 2.5.5 build targets
#                                  Python 3.8 and will be silently skipped if
#                                  python is pinned to 3.10. Let the solver pick
#                                  the latest compatible build automatically.
#   - python                     : NOT pinned — let solver choose whatever
#                                  satisfies orthofinder + busco together.
# =============================================================================

cat > "${SCRIPT_DIR}/config/environment.yml" << 'YMLEOF'
name: cestoda_busco_builder

channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python>=3.8
  - busco=6.0.0
  - sepp=4.5.5
  - miniprot
  - metaeuk
  - hmmer>=3.3.2
  - orthofinder
  - diamond>=2.1
  - mafft>=7.520
  - fasttree>=2.1
  - biopython>=1.80
  - pandas>=1.5
  - numpy>=1.23
  - wget
  - curl
  - ncbi-datasets-cli
  - entrez-direct
  - pigz
  - seqkit
YMLEOF

echo "[OK] config/environment.yml written."

# =============================================================================
# STEP 3 — Source the config we just wrote
# =============================================================================

source "${SCRIPT_DIR}/config/project.cfg"

LOG="${LOG_DIR}/00_setup.log"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  Cestoda BUSCO lineage builder — environment setup"
echo "  $(date)"
echo "============================================================"
echo ""

# =============================================================================
# STEP 4 — Create full project directory tree
# =============================================================================

echo "[INFO] Creating project directory structure..."

for D in \
    "${RAW_PROTEOMES}" \
    "${CLEAN_PROTEOMES}" \
    "${FILTERED_PROTEOMES}" \
    "${OG_DIR}" \
    "${MARKER_DIR}" \
    "${ALIGN_DIR}" \
    "${HMM_DIR}" \
    "${HMMSEARCH_DIR}" \
    "${ODB_DIR}/hmms" \
    "${ODB_DIR}/info" \
    "${VALIDATION_DIR}" \
    "${LOG_DIR}" \
    "${SCRIPTS_DIR}"
do
    mkdir -p "${D}"
done

echo "[OK]   Directory tree created under: ${PROJECT_ROOT}/"

# =============================================================================
# STEP 5 — Verify conda is available
# =============================================================================

echo ""
echo "[INFO] Checking conda..."

if ! command -v conda &>/dev/null; then
    echo "[ERROR] conda not found in PATH."
    echo "        Make sure Anaconda is installed and your shell is initialised."
    echo "        Try: source ~/anaconda3/etc/profile.d/conda.sh"
    exit 1
fi

echo "[OK]   $(conda --version)"

# =============================================================================
# STEP 6 — Create the conda environment
# =============================================================================

ENV_YML="${SCRIPT_DIR}/config/environment.yml"

echo ""
echo "[INFO] Creating environment '${CONDA_ENV}' from config/environment.yml..."
echo "       This may take 5–15 minutes."
echo ""

if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}"; then
    echo "       Environment '${CONDA_ENV}' already exists."
    read -rp "       Remove and recreate from scratch? [y/N] " REPLY
    if [[ "${REPLY}" =~ ^[Yy]$ ]]; then
        conda env remove -n "${CONDA_ENV}" --yes
        conda env create -f "${ENV_YML}"
        echo "[OK]   Environment recreated."
    else
        echo "       Skipping — using existing environment as-is."
    fi
else
    conda env create -f "${ENV_YML}"
    echo "[OK]   Environment '${CONDA_ENV}' created."
fi

# =============================================================================
# STEP 7 — Verify all required tools are functional
# =============================================================================

echo ""
echo "[INFO] Verifying installed tools..."
echo ""

# LD_LIBRARY_PATH must be exported before activation to avoid a bug in the
# aster package's conda activate script, which unconditionally appends to
# LD_LIBRARY_PATH without guarding against the variable being unset.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"

declare -A TOOL_TESTS=(
    ["busco"]="busco --version"
    ["datasets"]="datasets --version"
    ["diamond"]="diamond version"
    ["fasttree"]="FastTree -help"
    ["hmmbuild"]="hmmbuild -h"
    ["hmmemit"]="hmmemit -h"
    ["hmmsearch"]="hmmsearch -h"
    ["mafft"]="mafft --version"
    ["miniprot"]="miniprot --version"
    ["orthofinder"]="orthofinder --version"
    ["pigz"]="pigz --version"
    ["seqkit"]="seqkit version"
)

FAIL_COUNT=0
for TOOL in $(echo "${!TOOL_TESTS[@]}" | tr ' ' '\n' | sort); do
    if ${TOOL_TESTS[$TOOL]} &>/dev/null 2>&1; then
        printf "       %-20s  ✓\n" "${TOOL}"
    else
        printf "       %-20s  ✗  [NOT FOUND]\n" "${TOOL}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

echo ""
if [[ "${FAIL_COUNT}" -gt 0 ]]; then
    echo "[WARN] ${FAIL_COUNT} tool(s) failed."
    echo "       Remove and recreate the environment:"
    echo "         conda env remove -n ${CONDA_ENV} --yes"
    echo "         bash 00_setup_environment.sh"
else
    echo "[OK]   All tools verified."
fi

# =============================================================================
# Done
# =============================================================================

echo ""
echo "============================================================"
echo "  Setup complete."
echo ""
echo "  Activate the environment:"
echo "    conda activate ${CONDA_ENV}"
echo ""
echo "  Run the pipeline in order:"
echo "    bash 01_retrieve_data.sh"
echo "    bash 02_run_analysis.sh    # 8–24 h"
echo "    bash 03_build_odb.sh"
echo ""
echo "  CPU threads currently set to: ${OF_THREADS}"
echo "  To change: nano config/project.cfg  →  edit OF_THREADS"
echo ""
echo "  Log: ${LOG}"
echo "============================================================"
echo ""
