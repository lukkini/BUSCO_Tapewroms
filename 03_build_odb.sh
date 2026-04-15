#!/usr/bin/env bash
# =============================================================================
# 03_build_odb.sh
# -----------------------------------------------------------------------------
# Assembles the final BUSCO lineage dataset (ODB) directory and validates it.
#
#   Stage 1 — Ancestral sequence generation (hmmemit)
#   Stage 2 — Score cutoff computation (hmmsearch self-search)
#   Stage 3 — Compressed reference protein database (refseq_db.faa.gz)
#   Stage 4 — Dataset metadata files (dataset.cfg, info/, links)
#   Stage 5 — Structural integrity verification
#   Stage 6 — Self-validation on training Cestoda proteomes
#   Stage 7 — (Optional) Evaluation on a new genome assembly
#
# Output:
#   cestoda_odb_custom/   — ready for use with: busco -l cestoda_odb_custom/
#
# Usage (with environment activated):
#   bash 03_build_odb.sh
#
# To evaluate a new genome, set NEW_GENOME below or pass it as an argument:
#   bash 03_build_odb.sh /path/to/new_genome.fa
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"

# Guard against aster's conda activate script failing on unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Optional: accept a new genome path as a command-line argument
NEW_GENOME="${NEW_GENOME:-${1:-}}"

mkdir -p "${HMMSEARCH_DIR}" "${ODB_DIR}/hmms" "${ODB_DIR}/info" \
         "${VALIDATION_DIR}" "${LOG_DIR}" "${SCRIPTS_DIR}"

LOG="${LOG_DIR}/03_build_odb.log"
exec > >(tee -a "${LOG}") 2>&1

# =============================================================================
# Runtime resource detection
# =============================================================================
# The config file provides default parallelism values, but the script adapts to
# the actual machine where it is being executed. Effective runtime values are
# capped so the script never oversubscribes available CPU resources.
if command -v nproc &>/dev/null; then
    AVAILABLE_CORES="$(nproc)"
else
    AVAILABLE_CORES="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi

if [[ -z "${AVAILABLE_CORES}" ]] || ! [[ "${AVAILABLE_CORES}" =~ ^[0-9]+$ ]] || [[ "${AVAILABLE_CORES}" -lt 1 ]]; then
    AVAILABLE_CORES=1
fi

EFFECTIVE_THREADS_PER_JOB="${THREADS_PER_JOB}"
if [[ "${EFFECTIVE_THREADS_PER_JOB}" -gt "${AVAILABLE_CORES}" ]]; then
    EFFECTIVE_THREADS_PER_JOB="${AVAILABLE_CORES}"
fi
if [[ "${EFFECTIVE_THREADS_PER_JOB}" -lt 1 ]]; then
    EFFECTIVE_THREADS_PER_JOB=1
fi

EFFECTIVE_PARALLEL_JOBS="${PARALLEL_JOBS}"
MAX_JOBS_BY_CORES=$(( AVAILABLE_CORES / EFFECTIVE_THREADS_PER_JOB ))
if [[ "${MAX_JOBS_BY_CORES}" -lt 1 ]]; then
    MAX_JOBS_BY_CORES=1
fi
if [[ "${EFFECTIVE_PARALLEL_JOBS}" -gt "${MAX_JOBS_BY_CORES}" ]]; then
    EFFECTIVE_PARALLEL_JOBS="${MAX_JOBS_BY_CORES}"
fi
if [[ "${EFFECTIVE_PARALLEL_JOBS}" -lt 1 ]]; then
    EFFECTIVE_PARALLEL_JOBS=1
fi

EFFECTIVE_PIGZ_THREADS="${OF_THREADS}"
if [[ "${EFFECTIVE_PIGZ_THREADS}" -gt "${AVAILABLE_CORES}" ]]; then
    EFFECTIVE_PIGZ_THREADS="${AVAILABLE_CORES}"
fi
if [[ "${EFFECTIVE_PIGZ_THREADS}" -lt 1 ]]; then
    EFFECTIVE_PIGZ_THREADS=1
fi

echo ""
echo "============================================================"
echo "  03_build_odb.sh — $(date)"
echo "============================================================"
echo ""
echo "[INFO] Runtime resource detection"
echo "  Available CPU cores          : ${AVAILABLE_CORES}"
echo "  Requested PARALLEL_JOBS      : ${PARALLEL_JOBS}"
echo "  Requested THREADS_PER_JOB    : ${THREADS_PER_JOB}"
echo "  Effective parallel jobs      : ${EFFECTIVE_PARALLEL_JOBS}"
echo "  Effective threads per job    : ${EFFECTIVE_THREADS_PER_JOB}"
echo "  Effective pigz threads       : ${EFFECTIVE_PIGZ_THREADS}"
echo ""

# ── Pre-flight checks ─────────────────────────────────────────────────────────
N_HMM=$(find "${HMM_DIR}" -name "*.hmm" -size +0c 2>/dev/null | wc -l)
if [[ "${N_HMM}" -eq 0 ]]; then
    echo "[ERROR] No HMM profiles found in ${HMM_DIR}/"
    echo "        Run 02_run_analysis.sh first."
    exit 1
fi
echo "[INFO] ${N_HMM} HMM profiles found."

# =============================================================================
# Write the score-cutoff Python helper to disk
# =============================================================================

cat > "${SCRIPTS_DIR}/compute_score_cutoffs.py" << 'PYEOF'
#!/usr/bin/env python3
"""
compute_score_cutoffs.py
Derive per-marker hmmsearch bit-score cutoffs from self-hit results.

Method:
    For each marker profile, collect all full-sequence bit scores from an
    hmmsearch of that profile against the complete reference protein set.
    The cutoff is set to (minimum bit score) × FRACTION.

    This means BUSCO will accept as a true hit any sequence scoring ≥ FRACTION
    of the worst-scoring confirmed true ortholog in the training set, giving a
    FRACTION-sized buffer to accommodate divergent copies in new genomes.

Usage:
    python compute_score_cutoffs.py \
        --hmmsearch_dir <dir> \
        --output        <scores_cutoff_path> \
        --fraction      <float, default 0.90>
"""

import argparse
import glob
import os
import sys


def parse_tblout(path: str) -> list[float]:
    """
    Parse an hmmsearch --tblout file.
    Returns full-sequence bit scores for all significant hits.
    Lines beginning with '#' are comments; data columns are whitespace-separated.
    Column indices (0-based): 0=target, 2=query, 4=E-value, 5=score.
    """
    scores = []
    with open(path) as fh:
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
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hmmsearch_dir", required=True,
                   help="Directory containing .tbl files from hmmsearch")
    p.add_argument("--output",        required=True,
                   help="Output path for scores_cutoff file")
    p.add_argument("--fraction",      type=float, default=0.90,
                   help="Multiplier applied to minimum hit score (default: 0.90)")
    args = p.parse_args()

    tbl_files = glob.glob(os.path.join(args.hmmsearch_dir, "*.tbl"))
    if not tbl_files:
        print(f"[ERROR] No .tbl files found in {args.hmmsearch_dir}")
        sys.exit(1)

    cutoffs: dict[str, float] = {}
    no_hits: list[str] = []

    for tbl in sorted(tbl_files):
        og = os.path.basename(tbl).replace(".tbl", "")
        scores = parse_tblout(tbl)

        if not scores:
            # Profile matched nothing — assign an effectively unreachable cutoff.
            # These markers will appear "missing" in all BUSCO runs, which is
            # the correct behaviour: we cannot reliably use a profile with zero
            # reference hits.
            no_hits.append(og)
            cutoffs[og] = 999.0
        else:
            cutoffs[og] = round(min(scores) * args.fraction, 2)

    with open(args.output, "w") as out:
        for og in sorted(cutoffs):
            out.write(f"{og}\t{cutoffs[og]}\n")

    print(f"  Written: {len(cutoffs)} cutoffs → {args.output}")

    if no_hits:
        print(f"\n  WARNING: {len(no_hits)} profiles produced no self-hits.")
        print("    These markers will always appear 'missing' in BUSCO runs.")
        print("    Consider removing them or re-examining the alignment quality.")
        for og in no_hits[:5]:
            print(f"    - {og}")
        if len(no_hits) > 5:
            print(f"    ... and {len(no_hits) - 5} more (see log for full list)")


if __name__ == "__main__":
    main()
PYEOF

chmod +x "${SCRIPTS_DIR}/compute_score_cutoffs.py"

# =============================================================================
# STAGE 1 — Ancestral sequence generation (hmmemit)
# =============================================================================
# BUSCO stores two auxiliary sequence files used by BLAST-based genome search:
#
#   ancestral
#     One consensus sequence per marker, emitted with `hmmemit -c`.
#     The consensus sequence is the most probable amino acid at each match
#     state of the HMM profile — equivalent to the "consensus string" of the
#     alignment used to build the profile.  Used for initial tBLASTn searches
#     in Augustus and GeneMark pipelines.
#
#   ancestral_variants
#     Five sampled sequences per marker, emitted with `hmmemit -N 5`.
#     These stochastically sample the profile distribution, capturing
#     sequence diversity around the consensus.  Used to improve sensitivity
#     when the consensus is too divergent from the target genome's sequences.
#
# Both files are FASTA-formatted with one header per record.
# Headers are normalised to the orthogroup ID (BUSCO requirement).
#
# Note: Miniprot and Metaeuk (BUSCO v6 default) do NOT use these files.
#       They are generated for compatibility with all BUSCO run modes.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Ancestral sequence generation (hmmemit)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

ANCESTRAL="${ODB_DIR}/ancestral"
ANCESTRAL_VARS="${ODB_DIR}/ancestral_variants"

# Always regenerate to ensure consistency with current HMM set
: > "${ANCESTRAL}"
: > "${ANCESTRAL_VARS}"

N_EMITTED=0
for HMM in "${HMM_DIR}"/*.hmm; do
    OG=$(basename "${HMM}" .hmm)

    # Consensus sequence: rename header to bare orthogroup ID
    hmmemit -c "${HMM}" 2>/dev/null \
        | awk -v id="${OG}" '/^>/{print ">"id; next} {print}' \
        >> "${ANCESTRAL}"

    # Variant sequences: rename headers to <OG>_1 .. <OG>_5
    hmmemit -N 5 "${HMM}" 2>/dev/null \
        | awk -v id="${OG}" '
            BEGIN { n = 0 }
            /^>/  { n++; print ">"id"_"n; next }
                  { print }
        ' >> "${ANCESTRAL_VARS}"

    N_EMITTED=$((N_EMITTED + 1))
done

N_CONS=$(grep -c '^>' "${ANCESTRAL}")
N_VARS=$(grep -c '^>' "${ANCESTRAL_VARS}")
echo "  Profiles processed       : ${N_EMITTED}"
echo "  ancestral sequences      : ${N_CONS}"
echo "  ancestral_variants seqs  : ${N_VARS}"

# =============================================================================
# STAGE 2 — Score cutoff computation
# =============================================================================
# Run each HMM profile against ALL reference proteins (Cestoda + outgroups).
# Using the combined reference set (rather than Cestoda only) means the cutoff
# reflects the minimum score needed to detect the most divergent true ortholog
# in any training species, making it robust across the training-set range.
#
# The self-search uses E-value ≤ 1e-5 to restrict to biologically plausible
# hits while excluding spurious low-complexity matches.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 2 — Score cutoff computation (hmmsearch)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

ALL_PROTEINS="${HMMSEARCH_DIR}/all_reference_proteins.faa"

if [[ ! -s "${ALL_PROTEINS}" ]]; then
    echo "  Concatenating all reference proteomes..."
    cat "${FILTERED_PROTEOMES}"/*.fa > "${ALL_PROTEINS}"
    N_REFSEQ=$(grep -c '^>' "${ALL_PROTEINS}")
    echo "  Reference DB: ${N_REFSEQ} proteins"
else
    N_REFSEQ=$(grep -c '^>' "${ALL_PROTEINS}")
    echo "  Reference DB already exists: ${N_REFSEQ} proteins"
fi

# Build a DIAMOND database for faster repeated searches (optional; hmmsearch
# does not use DIAMOND — this is for reference only)

# Run hmmsearch in a parallel job pool
echo ""
echo "  Running hmmsearch: ${EFFECTIVE_PARALLEL_JOBS} parallel jobs × ${EFFECTIVE_THREADS_PER_JOB} threads"
echo ""

run_hmmsearch_pool() {
    local PIDS=()
    local N_SEARCHED=0
    local N_SKIP=0
    local N_FAIL=0

    for HMM in "${HMM_DIR}"/*.hmm; do
        OG=$(basename "${HMM}" .hmm)
        TBL="${HMMSEARCH_DIR}/${OG}.tbl"

        if [[ -f "${TBL}" ]] && [[ -s "${TBL}" ]]; then
            N_SKIP=$((N_SKIP + 1))
            continue
        fi

        hmmsearch \
            --cpu "${EFFECTIVE_THREADS_PER_JOB}" \
            --tblout "${TBL}" \
            -E 1e-5 \
            --noali \
            "${HMM}" \
            "${ALL_PROTEINS}" \
            > /dev/null 2>>"${LOG_DIR}/hmmsearch.log" &

        PIDS+=($!)
        N_SEARCHED=$((N_SEARCHED + 1))

        if [[ "${#PIDS[@]}" -ge "${EFFECTIVE_PARALLEL_JOBS}" ]]; then
            for PID in "${PIDS[@]}"; do
                wait "${PID}" || N_FAIL=$((N_FAIL + 1))
            done
            PIDS=()
            echo "  Searched ${N_SEARCHED} profiles so far..."
        fi
    done

    for PID in "${PIDS[@]}"; do
        wait "${PID}" || N_FAIL=$((N_FAIL + 1))
    done

    echo "  hmmsearch: ${N_SEARCHED} new / ${N_SKIP} already done / ${N_FAIL} failed"

    if [[ "${N_FAIL}" -ne 0 ]]; then
        echo "  [ERROR] ${N_FAIL} hmmsearch job(s) failed — check ${LOG_DIR}/hmmsearch.log"
        return 1
    fi
}

run_hmmsearch_pool

N_TBL=$(find "${HMMSEARCH_DIR}" -name "*.tbl" -size +0c | wc -l)
echo "  hmmsearch tables present: ${N_TBL}"

if [[ "${N_TBL}" -ne "${N_HMM}" ]]; then
    echo "  [ERROR] Expected ${N_HMM} hmmsearch table(s), found ${N_TBL}"
    echo "          Check ${LOG_DIR}/hmmsearch.log and ${HMMSEARCH_DIR}/"
    exit 1
fi

echo ""
echo "  Computing score cutoffs (fraction: ${SCORE_CUTOFF_FRACTION})..."
python "${SCRIPTS_DIR}/compute_score_cutoffs.py" \
    --hmmsearch_dir "${HMMSEARCH_DIR}" \
    --output        "${ODB_DIR}/scores_cutoff" \
    --fraction      "${SCORE_CUTOFF_FRACTION}"

N_CUTOFFS=$(wc -l < "${ODB_DIR}/scores_cutoff")
echo "  scores_cutoff: ${N_CUTOFFS} entries"

if [[ "${N_CUTOFFS}" -ne "${N_HMM}" ]]; then
    echo "  [ERROR] Expected ${N_HMM} score cutoff entries, found ${N_CUTOFFS}"
    echo "          scores_cutoff must contain exactly one row per HMM."
    exit 1
fi

# =============================================================================
# STAGE 3 — Reference protein database
# =============================================================================
# refseq_db.faa.gz is the protein database searched by Miniprot and Metaeuk
# when BUSCO runs in genome mode.  These tools map reference proteins onto
# the new assembly to generate candidate gene models for HMM scoring.
#
# IMPORTANT: We include only Cestoda proteins here.
#   Including outgroup proteins (Trematoda, Monogenea, etc.) would cause
#   BUSCO to map non-Cestoda proteins onto the new genome, potentially
#   generating spurious gene models that partially satisfy HMM profiles
#   and inflate apparent completeness scores.
#
# The database is compressed with pigz (parallel gzip) if available,
# otherwise falls back to gzip.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 3 — Reference protein database"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

REFSEQ_DB="${ODB_DIR}/refseq_db.faa.gz"

if [[ -f "${REFSEQ_DB}" ]] && [[ -s "${REFSEQ_DB}" ]]; then
    echo "  [SKIP] refseq_db.faa.gz already exists."
else
    TMPDB="${ODB_DIR}/refseq_db.faa"
    : > "${TMPDB}"

    echo "  Adding Cestoda reference proteins:"
    for TAG in "${CESTODA_TAGS[@]}"; do
        FA="${FILTERED_PROTEOMES}/${TAG}.fa"
        if [[ -f "${FA}" ]] && [[ -s "${FA}" ]]; then
            N=$(grep -c '^>' "${FA}")
            cat "${FA}" >> "${TMPDB}"
            printf "    %-48s  %d proteins\n" "${TAG}" "${N}"
        else
            echo "    [WARN] Missing: ${TAG}.fa — skipping"
        fi
    done

    N_REFDB=$(grep -c '^>' "${TMPDB}")
    echo ""
    echo "  Total: ${N_REFDB} proteins"
    echo "  Compressing..."

    if command -v pigz &>/dev/null; then
        pigz -9 -p "${EFFECTIVE_PIGZ_THREADS}" "${TMPDB}"
    else
        gzip -9 "${TMPDB}"
    fi

    echo "  refseq_db.faa.gz: $(du -sh "${REFSEQ_DB}" | cut -f1)"
fi

# Copy HMM profiles into the ODB directory
echo ""
echo "  Syncing HMM profiles to ODB..."
mkdir -p "${ODB_DIR}/hmms"
if command -v rsync &>/dev/null; then
    rsync -a --include="*.hmm" --exclude="*" "${HMM_DIR}/" "${ODB_DIR}/hmms/"
else
    find "${ODB_DIR}/hmms" -name "*.hmm" -delete
    cp "${HMM_DIR}"/*.hmm "${ODB_DIR}/hmms/"
fi
echo "  Profiles in ODB: $(find "${ODB_DIR}/hmms" -name "*.hmm" | wc -l)"

# =============================================================================
# STAGE 4 — Dataset metadata
# =============================================================================
# BUSCO reads the following files at runtime:
#
#   dataset.cfg
#     Key-value configuration that identifies the lineage, specifies the
#     number of markers, and sets the biological domain (eukaryota triggers
#     the Miniprot/Metaeuk/Augustus gene-prediction pipeline).
#
#   info/ogs.id.info
#     Newline-separated list of all BUSCO marker IDs (orthogroup IDs).
#     BUSCO uses this to know which genes to look for and to compute the
#     "Missing" fraction = (listed genes not found) / (total listed genes).
#
#   info/species.info
#     Documentation of which species were used to build the dataset.
#     Not used programmatically by BUSCO but required for provenance.
#
#   links_to_ODB.txt
#     Annotation file linking marker IDs to functional information.
#     In official BUSCO datasets this links to OrthoDB identifiers.
#     We populate it with our conservation-tag metadata instead.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 4 — Dataset metadata"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

N_BUSCOS=$(find "${ODB_DIR}/hmms" -name "*.hmm" | wc -l)
TODAY=$(date +%Y-%m-%d)

# ── dataset.cfg ───────────────────────────────────────────────────────────────
cat > "${ODB_DIR}/dataset.cfg" << EOF
[busco_dataset_info]
name = ${DATASET_NAME}
lineage = ${DATASET_LINEAGE}
creation_date = ${TODAY}
number_of_buscos = ${N_BUSCOS}
number_of_species = 11
domain = ${DATASET_DOMAIN}
EOF

echo "  Written: dataset.cfg"

# ── info/ogs.id.info ──────────────────────────────────────────────────────────
find "${ODB_DIR}/hmms" -name "*.hmm" \
    | xargs -I{} basename {} .hmm \
    | sort > "${ODB_DIR}/info/ogs.id.info"
echo "  Written: info/ogs.id.info  (${N_BUSCOS} entries)"

# ── info/species.info ─────────────────────────────────────────────────────────
cat > "${ODB_DIR}/info/species.info" << 'EOF'
# Custom Cestoda BUSCO dataset — species provenance
# species_tag  clade  bioproject  role  assembly_notes
echinococcus_canadensis_g7      Cestoda       PRJEB8992    target
echinococcus_granulosus_g1      Cestoda       PRJEB121     target
echinococcus_multilocularis     Cestoda       PRJEB122     target
hymenolepis_diminuta            Cestoda       PRJEB30942   target
hymenolepis_microstoma          Cestoda       PRJEB124     target
hymenolepis_nana                Cestoda       PRJEB508     target    fragmented_N50=19kb
mesocestoides_corti             Cestoda       PRJEB510     target    fragmented_N50=65kb
taenia_asiatica                 Cestoda       PRJEB532     target
taenia_multiceps                Cestoda       PRJNA307624  target
taenia_saginata                 Cestoda       PRJNA71493   target
taenia_solium                   Cestoda       PRJNA170813  target
clonorchis_sinensis             Trematoda     PRJNA386618  outgroup
fasciola_hepatica               Trematoda     PRJEB58756   outgroup
heterobilharzia_americana      Trematoda     TD1_PRJEB44434   outgroup
opisthorchis_felineus           Trematoda     PRJNA413383  outgroup
opisthorchis_viverrini          Trematoda     PRJNA222628  outgroup
paragonimus_westermani          Trematoda     PRJNA219632  outgroup
schistosoma_mansoni             Trematoda     PRJEA36577   outgroup
schistosoma_haematobium         Trematoda     PRJNA78265   outgroup
schistosoma_japonicum           Trematoda     PRJNA520774  outgroup
schistosoma_rodhaini            Trematoda     TD1_PRJEB44434   outgroup
trichobilharzia_regenti         Trematoda     PRJEB44434   outgroup
trichobilharzia_szidati         Trematoda     PRJEB44434   outgroup
schmidtea_mediterranea          Rhabditophora S2F19H1_PRJNA885486  outgroup  schMedS3_haplotype1
gyrodactylus_bullatarudis       Monogenea     PRJNA532341  outgroup
gyrodactylus_salaris            Monogenea     PRJNA244375  outgroup  fragmented_N50=18kb
EOF
echo "  Written: info/species.info"

# ── links_to_ODB.txt ─────────────────────────────────────────────────────────
META="${MARKER_DIR}/marker_metadata.tsv"
LINKS="${ODB_DIR}/links_to_ODB.txt"

echo "# Custom Cestoda BUSCO dataset — marker metadata" > "${LINKS}"
echo "# orthogroup_id  n_cestoda_present  n_cestoda_singlecopy  singlecopy_frac  n_outgroup_present  conservation" >> "${LINKS}"

if [[ -f "${META}" ]]; then
    tail -n +2 "${META}" >> "${LINKS}"
    echo "  Written: links_to_ODB.txt  ($(wc -l < "${LINKS}") lines)"
else
    echo "  [WARN] marker_metadata.tsv not found — links_to_ODB.txt will have headers only"
fi

# =============================================================================
# STAGE 5 — Structural integrity check
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 5 — Dataset integrity verification"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

REQUIRED=(
    "hmms"
    "ancestral"
    "ancestral_variants"
    "scores_cutoff"
    "refseq_db.faa.gz"
    "dataset.cfg"
    "info/ogs.id.info"
    "info/species.info"
    "links_to_ODB.txt"
)

FAIL=0
for ITEM in "${REQUIRED[@]}"; do
    FULL="${ODB_DIR}/${ITEM}"
    if [[ -e "${FULL}" ]] && [[ -s "${FULL}" ]]; then
        printf "  ✓  %-30s\n" "${ITEM}"
    else
        printf "  ✗  %-30s  [MISSING OR EMPTY]\n" "${ITEM}"
        FAIL=1
    fi
done

# Cross-check: every ID in ogs.id.info should have a corresponding .hmm file
echo ""
echo "  Cross-checking IDs vs HMM files..."
MISSING_HMM=0
while IFS= read -r OG_ID; do
    if [[ ! -f "${ODB_DIR}/hmms/${OG_ID}.hmm" ]]; then
        echo "  [WARN] ${OG_ID} is in ogs.id.info but has no matching .hmm file"
        MISSING_HMM=$((MISSING_HMM + 1))
    fi
done < "${ODB_DIR}/info/ogs.id.info"

if [[ "${MISSING_HMM}" -eq 0 ]]; then
    echo "  ✓  All OG IDs have matching HMM profiles"
else
    echo "  ✗  ${MISSING_HMM} OG IDs have no matching HMM profile"
    FAIL=1
fi

if [[ "${N_BUSCOS}" -ne "${N_HMM}" ]]; then
    echo "  [ERROR] HMM count mismatch: ODB has ${N_BUSCOS}, source HMM_DIR has ${N_HMM}"
    FAIL=1
fi

if [[ "$(wc -l < "${ODB_DIR}/info/ogs.id.info")" -ne "${N_BUSCOS}" ]]; then
    echo "  [ERROR] info/ogs.id.info entry count does not match HMM count"
    FAIL=1
fi

if [[ "$(wc -l < "${ODB_DIR}/scores_cutoff")" -ne "${N_BUSCOS}" ]]; then
    echo "  [ERROR] scores_cutoff entry count does not match HMM count"
    FAIL=1
fi

if [[ "${FAIL}" -ne 0 ]]; then
    echo ""
    echo "  [ERROR] Dataset is incomplete. Review the items marked ✗ above."
    exit 1
fi

echo ""
echo "  Dataset is complete and internally consistent."
echo ""
echo "  Summary:"
printf "    %-35s %d\n" "BUSCO markers:" "${N_BUSCOS}"
printf "    %-35s %d\n" "Score cutoffs:" "$(wc -l < "${ODB_DIR}/scores_cutoff")"
printf "    %-35s %d\n" "Ancestral sequences:" "${N_CONS}"
printf "    %-35s %s\n" "Reference DB size:" "$(du -sh "${ODB_DIR}/refseq_db.faa.gz" | cut -f1)"

# =============================================================================
# STAGE 6 — Self-validation on training Cestoda proteomes
# =============================================================================
# Running BUSCO on species that were used to build the dataset gives an upper
# bound on expected completeness scores.  If the training species score <90%
# on their own dataset, there is a problem with the HMM profiles, score
# cutoffs, or the refseq_db.
#
# We use protein mode (-m proteins) for speed — gene prediction is not needed
# when we already have the annotated proteome.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 6 — Self-validation (training species)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "  Validation species: ${VALIDATION_SPECIES[*]}"
echo ""

for SP in "${VALIDATION_SPECIES[@]}"; do
    FA="${FILTERED_PROTEOMES}/${SP}.fa"
    OUT_NAME="val_${SP}"
    OUT_DIR="${VALIDATION_DIR}/${OUT_NAME}"

    if [[ ! -f "${FA}" ]]; then
        echo "  [SKIP] ${SP} — proteome not found: ${FA}"
        continue
    fi

    if [[ -d "${OUT_DIR}" ]]; then
        echo "  [SKIP] ${SP} — validation run already exists: ${OUT_DIR}"
    else
        echo "  [RUN]  ${SP}..."
        busco \
            -i  "${FA}" \
            -m  proteins \
            -l  "${ODB_DIR}" \
            -o  "${OUT_NAME}" \
            --out_path "${VALIDATION_DIR}" \
            -c  "${OF_THREADS}" \
            --offline \
            --quiet \
            2>>"${LOG_DIR}/busco_validation.log" || {
                echo "  [WARN] BUSCO exited non-zero for ${SP} — check ${LOG_DIR}/busco_validation.log"
            }
    fi

    # Print the short summary from the BUSCO output
    SUMMARY=$(find "${OUT_DIR}" -name "short_summary*.txt" 2>/dev/null | head -1)
    if [[ -f "${SUMMARY}" ]]; then
        echo ""
        echo "  ── ${SP} ──────────────────────────────────────────"
        grep -E "C:|F:|M:|Complete|Missing|Fragmented|Total" "${SUMMARY}" \
            | head -8 \
            | sed 's/^/    /'
        echo ""
    fi
done

# =============================================================================
# STAGE 7 — Optional: evaluation on a new genome assembly
# =============================================================================
# Usage examples:
#   bash 03_build_odb.sh /path/to/my_new_cestode.fa
#   NEW_GENOME=/path/to/my_new_cestode.fa bash 03_build_odb.sh
#
# Expected completeness for a good Cestoda genome: 65–85%.
# Tapeworms have undergone extensive gene loss (34+ homeobox families,
# multiple metabolic pathways) relative to free-living flatworms.
# Scores below 65% should prompt inspection of assembly statistics (N50,
# genome size, k-mer completeness).
# =============================================================================

if [[ -n "${NEW_GENOME}" ]]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Stage 7 — New genome evaluation"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    if [[ ! -f "${NEW_GENOME}" ]]; then
        echo "  [ERROR] File not found: ${NEW_GENOME}"
        exit 1
    fi

    GENOME_NAME=$(basename "${NEW_GENOME}" | sed 's/\.[^.]*$//')
    GENOME_OUT="eval_${GENOME_NAME}"

    echo "  Input: ${NEW_GENOME}"
    echo "  Output: ${VALIDATION_DIR}/${GENOME_OUT}"
    echo ""

    busco \
        -i  "${NEW_GENOME}" \
        -m  genome \
        -l  "${ODB_DIR}" \
        -o  "${GENOME_OUT}" \
        --out_path "${VALIDATION_DIR}" \
        -c  "${OF_THREADS}" \
        --offline \
        2>>"${LOG_DIR}/busco_new_genome.log"

    SUMMARY=$(find "${VALIDATION_DIR}/${GENOME_OUT}" -name "short_summary*.txt" 2>/dev/null | head -1)
    if [[ -f "${SUMMARY}" ]]; then
        echo ""
        echo "  ── Results for ${GENOME_NAME} ─────────────────────"
        cat "${SUMMARY}"
    fi
else
    echo ""
    echo "  [INFO] No new genome specified for evaluation."
    echo "         To test a genome assembly, run:"
    echo ""
    echo "           bash 03_build_odb.sh /path/to/genome.fa"
    echo ""
    echo "         Or use BUSCO directly:"
    echo ""
    echo "           busco -i <genome.fa> \\"
    echo "                 -m genome \\"
    echo "                 -l ${ODB_DIR} \\"
    echo "                 -o my_run \\"
    echo "                 -c ${OF_THREADS}"
fi

# =============================================================================
# FINAL SUMMARY
# =============================================================================

echo ""
echo "════════════════════════════════════════════════════════"
echo "  Custom Cestoda BUSCO dataset successfully built!"
echo "  $(date)"
echo "════════════════════════════════════════════════════════"
echo ""
echo "  Dataset location : ${ODB_DIR}/"
echo "  BUSCO markers    : ${N_BUSCOS}"
echo ""
echo "  Interpreting BUSCO scores:"
echo "    85–100%  Excellent — assembly completeness comparable to best training genomes"
echo "    70–85%   Good — typical high-quality Cestoda genome"
echo "    60–70%   Acceptable — consistent with known gene loss in Cestoda"
echo "    <60%     Investigate — likely assembly fragmentation"
echo ""
echo "  Conservation tags for each marker:"
echo "    ${MARKER_DIR}/marker_metadata.tsv"
echo ""
echo "    platyhelminthes_conserved → low score may reflect assembly quality"
echo "    cestoda_specific          → low score may reflect genuine gene loss"
echo ""
echo "  Logs: ${LOG_DIR}/"
echo ""
