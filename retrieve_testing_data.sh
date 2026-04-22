#!/usr/bin/env bash
# =============================================================================
# retrieve_testing_data.sh
# -----------------------------------------------------------------------------
# Download validation assets for the SAME species used in lineage construction,
# without modifying 01_retrieve_data.sh.
#
# Downloads BOTH:
#   - protein FASTA  -> data/01_raw_proteomes/<tag>.fa
#   - masked genome  -> data/00_validation_genomes/<tag>.fa.gz
#
# Notes:
#   - Existing files are reused and skipped safely.
#   - Protein retrieval mirrors the species/project tuples used for lineage
#     construction, but lives in this separate validation helper so the original
#     01_retrieve_data.sh remains untouched.
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

VALIDATION_GENOMES="${PROJECT_ROOT}/data/00_validation_genomes"

mkdir -p "${RAW_PROTEOMES}" "${VALIDATION_GENOMES}" "${LOG_DIR}"

LOG="${LOG_DIR}/retrieve_testing_data.log"
PROT_MANIFEST="${LOG_DIR}/testing_species_manifest.tsv"
GENOME_MANIFEST="${LOG_DIR}/testing_genome_manifest.tsv"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  retrieve_testing_data.sh — $(date)"
echo "  Source : WormBase ParaSite ${WBPS_RELEASE}"
echo "  FTP    : ${WBPS_BASE}"
echo "============================================================"
echo ""

printf "tag\tclade\trole\tspecies_dir\tftp_bioproj\turl\tn_sequences\tstatus\n" > "${PROT_MANIFEST}"
printf "tag\tclade\trole\tspecies_dir\tftp_bioproj\turl\tstatus\n" > "${GENOME_MANIFEST}"

download_proteome_wbps() {
    local SPECIES_DIR="$1"
    local FTP_BIOPROJ="$2"
    local TAG="$3"
    local CLADE="$4"
    local ROLE="$5"

    local OUT="${RAW_PROTEOMES}/${TAG}.fa"
    local TMP="${RAW_PROTEOMES}/.${TAG}.protein.fa.gz.part"
    local FAILED_SENTINEL="${RAW_PROTEOMES}/${TAG}.FAILED"
    local URL="${WBPS_BASE}/${SPECIES_DIR}/${FTP_BIOPROJ}/${SPECIES_DIR}.${FTP_BIOPROJ}.${WBPS_RELEASE}.protein.fa.gz"

    if [[ -s "${OUT}" ]]; then
        local N_EXISTING
        N_EXISTING=$(grep -c '^>' "${OUT}" || true)
        printf "  [SKIP][protein] %-30s %8d sequences\n" "${TAG}" "${N_EXISTING}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "${N_EXISTING}" "skipped" >> "${PROT_MANIFEST}"
        return 0
    fi

    rm -f "${TMP}" "${FAILED_SENTINEL}"

    printf "  [DL]  [protein] %s\n" "${TAG}"
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
        echo "  [FAIL][protein] ${TAG} — all download attempts failed"
        rm -f "${TMP}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_download" >> "${PROT_MANIFEST}"
        return 1
    fi

    if ! gzip -t "${TMP}" >/dev/null 2>&1; then
        echo "  [FAIL][protein] ${TAG} — downloaded file is not a valid gzip archive"
        rm -f "${TMP}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_not_gzip" >> "${PROT_MANIFEST}"
        return 1
    fi

    if ! gunzip -c "${TMP}" > "${OUT}"; then
        echo "  [FAIL][protein] ${TAG} — decompression failed"
        rm -f "${TMP}" "${OUT}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_gunzip" >> "${PROT_MANIFEST}"
        return 1
    fi

    rm -f "${TMP}" "${FAILED_SENTINEL}"

    local N_SEQ
    N_SEQ=$(grep -c '^>' "${OUT}" || true)
    if [[ "${N_SEQ}" -eq 0 ]]; then
        echo "  [FAIL][protein] ${TAG} — decompressed FASTA contains zero sequences"
        rm -f "${OUT}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_empty_fasta" >> "${PROT_MANIFEST}"
        return 1
    fi

    printf "  [OK]  [protein] %-30s %8d sequences\n" "${TAG}" "${N_SEQ}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
        "${URL}" "${N_SEQ}" "ok" >> "${PROT_MANIFEST}"
    return 0
}

download_genome_wbps() {
    local SPECIES_DIR="$1"
    local FTP_BIOPROJ="$2"
    local TAG="$3"
    local CLADE="$4"
    local ROLE="$5"

    local OUT="${VALIDATION_GENOMES}/${TAG}.fa.gz"
    local TMP="${VALIDATION_GENOMES}/.${TAG}.genomic_masked.fa.gz.part"
    local FAILED_SENTINEL="${VALIDATION_GENOMES}/${TAG}.FAILED"
    local URL="${WBPS_BASE}/${SPECIES_DIR}/${FTP_BIOPROJ}/${SPECIES_DIR}.${FTP_BIOPROJ}.${WBPS_RELEASE}.genomic_masked.fa.gz"

    if [[ -s "${OUT}" ]]; then
        printf "  [SKIP][genome ] %-30s %s\n" "${TAG}" "$(basename "${OUT}")"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "skipped" >> "${GENOME_MANIFEST}"
        return 0
    fi

    rm -f "${TMP}" "${FAILED_SENTINEL}"

    printf "  [DL]  [genome ] %s\n" "${TAG}"
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
        echo "  [FAIL][genome ] ${TAG} — all download attempts failed"
        rm -f "${TMP}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_download" >> "${GENOME_MANIFEST}"
        return 1
    fi

    if ! gzip -t "${TMP}" >/dev/null 2>&1; then
        echo "  [FAIL][genome ] ${TAG} — downloaded file is not a valid gzip archive"
        rm -f "${TMP}"
        touch "${FAILED_SENTINEL}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
            "${URL}" "FAILED_not_gzip" >> "${GENOME_MANIFEST}"
        return 1
    fi

    mv "${TMP}" "${OUT}"
    rm -f "${FAILED_SENTINEL}"

    printf "  [OK]  [genome ] %-30s %s\n" "${TAG}" "$(basename "${OUT}")"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${TAG}" "${CLADE}" "${ROLE}" "${SPECIES_DIR}" "${FTP_BIOPROJ}" \
        "${URL}" "ok" >> "${GENOME_MANIFEST}"
    return 0
}

download_species_assets() {
    local SPECIES_DIR="$1"
    local FTP_BIOPROJ="$2"
    local TAG="$3"
    local CLADE="$4"
    local ROLE="$5"

    local STATUS=0
    download_proteome_wbps "${SPECIES_DIR}" "${FTP_BIOPROJ}" "${TAG}" "${CLADE}" "${ROLE}" || STATUS=1
    download_genome_wbps   "${SPECIES_DIR}" "${FTP_BIOPROJ}" "${TAG}" "${CLADE}" "${ROLE}" || STATUS=1
    return "${STATUS}"
}

FAILURES=0

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
    if ! download_species_assets ${SPEC}; then
        FAILURES=$((FAILURES + 1))
    fi
done

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
    if ! download_species_assets ${SPEC}; then
        FAILURES=$((FAILURES + 1))
    fi
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  C. Rhabditophora — 1 species (deep outgroup)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if ! download_species_assets \
    "schmidtea_mediterranea" \
    "S2F19H1_PRJNA885486" \
    "schmidtea_mediterranea" \
    "Rhabditophora" \
    "outgroup"; then
    FAILURES=$((FAILURES + 1))
fi

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
    if ! download_species_assets ${SPEC}; then
        FAILURES=$((FAILURES + 1))
    fi
done

echo ""
echo "============================================================"
echo "  Download summary"
echo "============================================================"
echo ""
echo "  Protein manifest : ${PROT_MANIFEST}"
echo "  Genome manifest  : ${GENOME_MANIFEST}"
echo "  Genome dir       : ${VALIDATION_GENOMES}"
echo ""

N_PROT_OK=$(awk -F'\t' 'NR>1 && $8=="ok"         {c++} END{print c+0}' "${PROT_MANIFEST}")
N_PROT_SKIP=$(awk -F'\t' 'NR>1 && $8=="skipped"  {c++} END{print c+0}' "${PROT_MANIFEST}")
N_PROT_FAIL=$(awk -F'\t' 'NR>1 && $8 ~ /^FAILED/ {c++} END{print c+0}' "${PROT_MANIFEST}")

N_GEN_OK=$(awk -F'\t' 'NR>1 && $7=="ok"         {c++} END{print c+0}' "${GENOME_MANIFEST}")
N_GEN_SKIP=$(awk -F'\t' 'NR>1 && $7=="skipped"  {c++} END{print c+0}' "${GENOME_MANIFEST}")
N_GEN_FAIL=$(awk -F'\t' 'NR>1 && $7 ~ /^FAILED/ {c++} END{print c+0}' "${GENOME_MANIFEST}")

echo "  Proteins:"
echo "    Successful downloads    : ${N_PROT_OK}"
echo "    Skipped existing files  : ${N_PROT_SKIP}"
echo "    Failures                : ${N_PROT_FAIL}"
echo ""
echo "  Genomes:"
echo "    Successful downloads    : ${N_GEN_OK}"
echo "    Skipped existing files  : ${N_GEN_SKIP}"
echo "    Failures                : ${N_GEN_FAIL}"
echo ""

if [[ "${FAILURES}" -gt 0 || "${N_PROT_FAIL}" -gt 0 || "${N_GEN_FAIL}" -gt 0 ]]; then
    echo "[WARN] One or more species assets failed. See manifests and *.FAILED sentinels."
    exit 1
fi

echo "[DONE] All requested proteins and genomes are available."
