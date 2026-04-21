#!/usr/bin/env bash
# =============================================================================
# 01_retrieve_data.sh
# -----------------------------------------------------------------------------
# Downloads protein FASTA files for all 26 species used to construct the
# custom Cestoda BUSCO lineage dataset from WormBase ParaSite WBPS19.
#
# Output:
#   data/01_raw_proteomes/<tag>.fa
#   logs/species_manifest.tsv
#   data/01_raw_proteomes/<tag>.FAILED    (sentinel on failure)
#
# Usage:
#   bash 01_retrieve_data.sh
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

mkdir -p "${RAW_PROTEOMES}" "${LOG_DIR}"
LOG="${LOG_DIR}/01_retrieve_data.log"
MANIFEST="${LOG_DIR}/species_manifest.tsv"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  01_retrieve_data.sh — $(date)"
echo "  Source : WormBase ParaSite ${WBPS_RELEASE}"
echo "  FTP    : ${WBPS_BASE}"
echo "============================================================"
echo ""

# Manifest columns:
#   tag, clade, role, species_dir, ftp_bioproj, url, n_sequences, status
printf "tag\tclade\trole\tspecies_dir\tftp_bioproj\turl\tn_sequences\tstatus\n" > "${MANIFEST}"

# -----------------------------------------------------------------------------
# download_wbps <species_dir> <ftp_bioproj> <tag> <clade> <role>
# -----------------------------------------------------------------------------
# species_dir : WormBase ParaSite species directory name
# ftp_bioproj : exact FTP path component (may be qualified, e.g. TD1_PRJEB44434)
# tag         : pipeline-wide unique species tag / output filename stem
# clade       : Cestoda / Trematoda / Monogenea / Rhabditophora
# role        : target / outgroup
# -----------------------------------------------------------------------------
download_wbps() {
    local SPECIES_DIR="$1"
    local FTP_BIOPROJ="$2"
    local TAG="$3"
    local CLADE="$4"
    local ROLE="$5"

    local OUT="${RAW_PROTEOMES}/${TAG}.fa"
    local TMP="${RAW_PROTEOMES}/.${TAG}.fa.gz.part"
    local FAILED_SENTINEL="${RAW_PROTEOMES}/${TAG}.FAILED"
    local URL="${WBPS_BASE}/${SPECIES_DIR}/${FTP_BIOPROJ}/${SPECIES_DIR}.${FTP_BIOPROJ}.${WBPS_RELEASE}.protein.fa.gz"

    # Skip already-downloaded proteomes.
    if [[ -s "${OUT}" ]]; then
        local N_EXISTING
        N_EXISTING=$(grep -c '^>' "${OUT}" || true)
        printf "  [SKIP] %-30s %8d sequences\n" "${TAG}" "${N_EXISTING}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "${N_EXISTING}" "skipped" >> "${MANIFEST}"
        return 0
    fi

    rm -f "${TMP}" "${FAILED_SENTINEL}"

    printf "  [DL]   %s\n" "${TAG}"
    printf "         URL: %s\n" "${URL}"

    local ATTEMPT=0
    local SUCCESS=0

    while [[ "${ATTEMPT}" -lt "${WGET_RETRIES}" ]]; do
        ATTEMPT=$((ATTEMPT + 1))

        if wget \
            --no-verbose \
            --show-progress \
            --retry-connrefused \
            --waitretry=10 \
            --tries=1 \
            --timeout="${WGET_TIMEOUT}" \
            -O "${TMP}" \
            "${URL}"; then
            SUCCESS=1
            break
        fi

        echo "         Attempt ${ATTEMPT}/${WGET_RETRIES} failed."
        sleep 5
    done

    if [[ "${SUCCESS}" -ne 1 ]]; then
        echo "  [FAIL] ${TAG} — all download attempts failed"
        echo "         Failed URL: ${URL}"
        rm -f "${TMP}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_download" >> "${MANIFEST}"
        return 1
    fi

    # Validate that the payload is a real gzip stream before decompressing.
    if ! gzip -t "${TMP}" >/dev/null 2>&1; then
        echo "  [FAIL] ${TAG} — downloaded file is not a valid gzip archive"
        echo "         Failed URL: ${URL}"
        rm -f "${TMP}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_not_gzip" >> "${MANIFEST}"
        return 1
    fi

    if ! gunzip -c "${TMP}" > "${OUT}"; then
        echo "  [FAIL] ${TAG} — decompression failed"
        echo "         Failed URL: ${URL}"
        rm -f "${TMP}" "${OUT}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_gunzip" >> "${MANIFEST}"
        return 1
    fi

    rm -f "${TMP}" "${FAILED_SENTINEL}"

    local N_SEQ
    N_SEQ=$(grep -c '^>' "${OUT}" || true)

    if [[ "${N_SEQ}" -eq 0 ]]; then
        echo "  [FAIL] ${TAG} — decompressed FASTA contains zero sequences"
        echo "         Failed URL: ${URL}"
        rm -f "${OUT}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_empty_fasta" >> "${MANIFEST}"
        return 1
    fi

    printf "  [OK]   %-30s %8d sequences\n" "${TAG}" "${N_SEQ}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
        "${URL}" "${N_SEQ}" "ok" >> "${MANIFEST}"

    return 0
}

FAILURES=0

# =============================================================================
# SECTION A — CESTODA (11 target species)
# =============================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  A. Cestoda — 11 species (target clade)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

for SPEC in \
    "echinococcus_canadensis PRJEB8992 echinococcus_canadensis_g7 Cestoda target" \
    "echinococcus_granulosus PRJEB121 echinococcus_granulosus_g1 Cestoda target" \
    "echinococcus_multilocularis PRJEB122 echinococcus_multilocularis Cestoda target" \
    "hymenolepis_diminuta PRJEB30942 hymenolepis_diminuta Cestoda target" \
    "hymenolepis_microstoma PRJEB124 hymenolepis_microstoma Cestoda target" \
    "hymenolepis_nana PRJEB508 hymenolepis_nana Cestoda target" \
    "mesocestoides_corti PRJEB510 mesocestoides_corti Cestoda target" \
    "taenia_asiatica PRJEB532 taenia_asiatica Cestoda target" \
    "taenia_multiceps PRJNA307624 taenia_multiceps Cestoda target" \
    "taenia_saginata PRJNA71493 taenia_saginata Cestoda target" \
    "taenia_solium PRJNA170813 taenia_solium Cestoda target"
do
    # shellcheck disable=SC2086
    if ! download_wbps ${SPEC}; then
        FAILURES=$((FAILURES + 1))
    fi
done

# =============================================================================
# SECTION B — TREMATODA (12 outgroup species)
# =============================================================================
# Verified FTP path resolution for previously failing species:
#   schistosoma_rodhaini    -> TD1_PRJEB44434
#   trichobilharzia_regenti -> PRJEB44434
#
# schistosoma_rodhaini has both TD1_PRJEB44434 and TD2_PRJEB44434 directories
# in WBPS19. TD1 is used here because it is the primary WBPS19 species page and
# exposes a valid protein FASTA.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  B. Trematoda — 12 species (outgroup)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

for SPEC in \
    "clonorchis_sinensis PRJNA386618 clonorchis_sinensis Trematoda outgroup" \
    "fasciola_hepatica PRJEB58756 fasciola_hepatica Trematoda outgroup" \
    "heterobilharzia_americana TD1_PRJEB44434 heterobilharzia_americana Trematoda outgroup" \
    "opisthorchis_felineus PRJNA413383 opisthorchis_felineus Trematoda outgroup" \
    "opisthorchis_viverrini PRJNA222628 opisthorchis_viverrini Trematoda outgroup" \
    "paragonimus_westermani PRJNA219632 paragonimus_westermani Trematoda outgroup" \
    "schistosoma_mansoni PRJEA36577 schistosoma_mansoni Trematoda outgroup" \
    "schistosoma_haematobium PRJNA78265 schistosoma_haematobium Trematoda outgroup" \
    "schistosoma_japonicum PRJNA520774 schistosoma_japonicum Trematoda outgroup" \
    "schistosoma_rodhaini TD1_PRJEB44434 schistosoma_rodhaini Trematoda outgroup" \
    "trichobilharzia_regenti PRJEB44434 trichobilharzia_regenti Trematoda outgroup" \
    "trichobilharzia_szidati PRJEB44434 trichobilharzia_szidati Trematoda outgroup"
do
    # shellcheck disable=SC2086
    if ! download_wbps ${SPEC}; then
        FAILURES=$((FAILURES + 1))
    fi
done

# =============================================================================
# SECTION C — RHABDITOPHORA (1 deep outgroup)
# =============================================================================
# schmidtea_mediterranea has three WBPS19 genome project directories:
#   PRJNA12585
#   S2F19H1_PRJNA885486
#   S2F19H2_PRJNA885486
# We use S2F19H1_PRJNA885486 here because it corresponds to
# schMedS3_haplotype1, which is the more complete assembly and therefore the
# preferred deep outgroup proteome for orthogroup inference.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  C. Rhabditophora — 1 species (deep outgroup)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if ! download_wbps \
    "schmidtea_mediterranea" \
    "S2F19H1_PRJNA885486" \
    "schmidtea_mediterranea" \
    "Rhabditophora" \
    "outgroup"; then
    FAILURES=$((FAILURES + 1))
fi

# =============================================================================
# SECTION D — MONOGENEA (2 outgroup species)
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  D. Monogenea — 2 species (outgroup)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

for SPEC in \
    "gyrodactylus_bullatarudis PRJNA532341 gyrodactylus_bullatarudis Monogenea outgroup" \
    "gyrodactylus_salaris PRJNA244375 gyrodactylus_salaris Monogenea outgroup"
do
    # shellcheck disable=SC2086
    if ! download_wbps ${SPEC}; then
        FAILURES=$((FAILURES + 1))
    fi
done

# =============================================================================
# Final summary
# =============================================================================

echo ""
echo "============================================================"
echo "  Download summary"
echo "============================================================"

aWK_BIN="$(command -v awk || true)"
if [[ -n "${aWK_BIN}" ]] && [[ -s "${MANIFEST}" ]]; then
    awk -F'\t' '
        NR==1 { next }
        {
            printf("  %-30s %8s  %s\n", $1, $7, $8)
        }
    ' "${MANIFEST}"
fi

N_OK=$(awk -F'\t' 'NR>1 && $8=="ok"      {c++} END{print c+0}' "${MANIFEST}")
N_SKIP=$(awk -F'\t' 'NR>1 && $8=="skipped" {c++} END{print c+0}' "${MANIFEST}")
N_FAIL=$(awk -F'\t' 'NR>1 && $8 ~ /^FAILED/ {c++} END{print c+0}' "${MANIFEST}")
N_TOTAL=$(awk -F'\t' 'NR>1 {c++} END{print c+0}' "${MANIFEST}")

echo ""
echo "  Total species processed : ${N_TOTAL}"
echo "  Successful downloads    : ${N_OK}"
echo "  Skipped existing files  : ${N_SKIP}"
echo "  Failures                : ${N_FAIL}"
echo "  Manifest                : ${MANIFEST}"
echo ""

if [[ "${FAILURES}" -gt 0 || "${N_FAIL}" -gt 0 ]]; then
    echo "[WARN] One or more species failed. See .FAILED sentinels and ${MANIFEST}."
    exit 1
fi

echo "[DONE] All requested species are available in ${RAW_PROTEOMES}/"
