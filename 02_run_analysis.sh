#!/usr/bin/env bash
# =============================================================================
# 02_run_analysis.sh
# -----------------------------------------------------------------------------
# Core analysis pipeline. Runs sequentially through five stages:
#
#   Stage 1 — Sequence quality control and cleaning
#             Remove artefactual sequences before orthogroup inference.
#
#   Stage 2 — Longest-isoform filtering
#             Keep one representative protein per gene locus.
#
#   Stage 3 — OrthoFinder
#             All-vs-all DIAMOND + MCL clustering to infer orthogroups.
#
#   Stage 4 — Marker selection
#             Clade-aware filtering to identify Cestoda BUSCO marker genes.
#
#   Stage 5 — Multiple sequence alignment (MAFFT)
#             Align each marker orthogroup for HMM profile construction.
#
#   Stage 6 — HMM profile construction (HMMER hmmbuild)
#             Build one profile HMM per marker.
#
# Estimated runtime: 8–24 hours (OrthoFinder dominates; scales with OF_THREADS).
#
# Usage (with environment activated):
#   bash 02_run_analysis.sh
#
# Resumability: each stage checks for existing output files and skips
# completed work. To restart a stage, delete its output directory.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/project.cfg"

# Guard against aster's conda activate script failing on unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}"

# Avoid literal "*.fa" / "*.msa" strings when a directory is empty.
shopt -s nullglob

mkdir -p "${CLEAN_PROTEOMES}" "${FILTERED_PROTEOMES}" \
         "${MARKER_DIR}" "${ALIGN_DIR}" "${HMM_DIR}" "${LOG_DIR}" "${SCRIPTS_DIR}"

LOG="${LOG_DIR}/02_analysis.log"
exec > >(tee -a "${LOG}") 2>&1

echo ""
echo "============================================================"
echo "  02_run_analysis.sh — $(date)"
echo "============================================================"
echo ""

# ── Verify raw proteomes exist before starting ────────────────────────────────
N_RAW=$(find "${RAW_PROTEOMES}" -name "*.fa" -size +0c 2>/dev/null | wc -l)
if [[ "${N_RAW}" -eq 0 ]]; then
    echo "[ERROR] No protein FASTA files found in ${RAW_PROTEOMES}/"
    echo "        Run 01_retrieve_data.sh first."
    exit 1
fi

if [[ "${N_RAW}" -lt 26 ]]; then
    echo "[WARN] Expected 26 proteomes but found ${N_RAW}."
    echo "       Continuing, but reduced taxon sampling may affect orthogroup inference."
fi

N_FAIL=$(find "${RAW_PROTEOMES}" -name "*.FAILED" 2>/dev/null | wc -l)
if [[ "${N_FAIL}" -gt 0 ]]; then
    echo "[WARN] ${N_FAIL} species failed to download (see .FAILED files)."
    echo "       Missing species may affect orthogroup detection."
    echo "       Continuing with ${N_RAW} available proteomes."
    echo ""
fi

# =============================================================================
# Write Python helper scripts to disk.
# All scripts use a quoted heredoc delimiter ('PYEOF') to prevent bash from
# expanding variables inside the Python source code.
# =============================================================================

echo "[INFO] Writing Python helper scripts..."

# -----------------------------------------------------------------------------
# clean_sequences.py
# -----------------------------------------------------------------------------
# Performs four cleaning operations on each proteome FASTA:
#
#   1. Remove sequences shorter than MIN_LENGTH amino acids.
#      Very short predictions (<30 aa) are almost always artefacts of
#      gene prediction errors (spurious ORFs, assembly gaps).
#
#   2. Remove sequences containing internal stop codons ('*' mid-sequence).
#      A '*' at the end of a sequence is a valid translation terminator;
#      a '*' in the middle indicates a prediction spanning a frameshift,
#      an assembly join, or a mis-annotated intron.
#
#   3. Remove sequences where the fraction of ambiguous residues ('X')
#      exceeds MAX_X_FRACTION.
#      High-X sequences arise from low-coverage assembly regions and are
#      unreliable for profile HMM construction.
#
#   4. Remove exact duplicate sequences (same MD5 of the sequence string).
#      Duplicates may arise from redundant protein entries or copy-paste
#      artefacts in genome annotation pipelines.
#
# A TSV report is written alongside the output FASTA documenting all filters.
# -----------------------------------------------------------------------------
cat > "${SCRIPTS_DIR}/clean_sequences.py" << 'PYEOF'
#!/usr/bin/env python3
"""
clean_sequences.py
Protein FASTA quality-control and cleaning.

Usage:
    python clean_sequences.py \
        --input  <in.fa> \
        --output <out.fa> \
        --report <qc_report.tsv> \
        --min_length  <int, default 30> \
        --max_x_frac  <float, default 0.10>
"""

import argparse
import hashlib
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def md5(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


def clean_fasta(
    input_path: str,
    output_path: str,
    report_path: str,
    min_length: int = 30,
    max_x_frac: float = 0.10,
) -> dict:
    """
    Clean a protein FASTA and return a dictionary of counts.
    """
    seen_hashes: set = set()
    kept: list[SeqRecord] = []

    counts = {
        "total":           0,
        "too_short":       0,
        "internal_stop":   0,
        "high_x":          0,
        "duplicate":       0,
        "kept":            0,
    }

    report_rows = []

    for rec in SeqIO.parse(input_path, "fasta"):
        counts["total"] += 1
        seq_str = str(rec.seq).upper()
        reason = None

        # Filter 1: minimum length
        if len(seq_str) < min_length:
            counts["too_short"] += 1
            reason = f"too_short ({len(seq_str)} aa)"

        # Filter 2: internal stop codons
        # Strip a trailing stop (valid translation terminator) before checking.
        elif "*" in seq_str.rstrip("*"):
            counts["internal_stop"] += 1
            pos = seq_str.index("*")
            reason = f"internal_stop (pos {pos})"

        # Filter 3: excessive ambiguous residues
        else:
            x_frac = seq_str.count("X") / len(seq_str)
            if x_frac > max_x_frac:
                counts["high_x"] += 1
                reason = f"high_X ({x_frac:.2%})"

        if reason is None:
            # Filter 4: exact duplicates (by sequence MD5)
            h = md5(seq_str.rstrip("*"))  # strip terminal stop before hashing
            if h in seen_hashes:
                counts["duplicate"] += 1
                reason = "duplicate"
            else:
                seen_hashes.add(h)
                # Remove trailing stop codon from the stored sequence
                clean_seq = seq_str.rstrip("*")
                rec.seq = rec.seq.__class__(clean_seq)
                kept.append(rec)
                counts["kept"] += 1

        report_rows.append((rec.id, len(seq_str), reason or "kept"))

    # Write output FASTA
    with open(output_path, "w") as fh:
        SeqIO.write(kept, fh, "fasta")

    # Write QC report
    with open(report_path, "w") as fh:
        fh.write("seq_id\tlength\tstatus\n")
        for row in report_rows:
            fh.write("\t".join(str(v) for v in row) + "\n")

    return counts


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",      required=True)
    p.add_argument("--output",     required=True)
    p.add_argument("--report",     required=True)
    p.add_argument("--min_length", type=int,   default=30)
    p.add_argument("--max_x_frac", type=float, default=0.10)
    args = p.parse_args()

    counts = clean_fasta(
        args.input, args.output, args.report,
        args.min_length, args.max_x_frac,
    )

    removed = counts["total"] - counts["kept"]
    pct_kept = 100 * counts["kept"] / max(counts["total"], 1)

    print(
        f"  {Path(args.input).name:<50}  "
        f"{counts['total']:>7} in  →  "
        f"{counts['kept']:>7} kept  ({pct_kept:.1f}%)  "
        f"[short:{counts['too_short']}  stop:{counts['internal_stop']}  "
        f"X:{counts['high_x']}  dup:{counts['duplicate']}]"
    )


if __name__ == "__main__":
    main()
PYEOF

# -----------------------------------------------------------------------------
# keep_longest_isoform.py
# -----------------------------------------------------------------------------
# Most genome annotations include multiple protein isoforms per gene locus.
# OrthoFinder's orthogroup inference operates at gene level, so including all
# isoforms creates two problems:
#   (a) A single gene appears as N copies → inflates apparent gene number.
#   (b) Different isoforms of the same gene may be split across orthogroups.
#
# We keep the longest protein per gene, which is the most common isoform in
# the reference databases and most likely to represent the full-length ORF.
#
# Isoform suffix patterns handled (most common across WormBase / NCBI / Ensembl):
#   WormBase:   gene.1, gene.2             → strip trailing .\d+
#   NCBI:       XP_001234.1                → strip trailing .\d+
#   Ensembl:    GENE-RA, GENE-RB           → strip trailing -R[A-Z]
#   GenBank:    gene-mRNA-1, gene-mRNA-2   → strip -mRNA-\d+
# -----------------------------------------------------------------------------
cat > "${SCRIPTS_DIR}/keep_longest_isoform.py" << 'PYEOF'
#!/usr/bin/env python3
"""
keep_longest_isoform.py
Retain only the longest protein sequence per gene locus.

Usage:
    python keep_longest_isoform.py <input.fa> <output.fa>
"""

import re
import sys
from Bio import SeqIO


ISOFORM_PATTERNS = [
    re.compile(r'\.\d+$'),                           # gene.1, XP_001234.1
    re.compile(r'[-_]m?[Rr][Nn][Aa][-_]?\d+$'),     # gene-mRNA-1, gene_rna_2
    re.compile(r'[-_][RT][A-Z]$'),                   # GENE-RA, GENE-RB (Ensembl)
    re.compile(r'[-_][Tt]ranscript[-_]?\d+$'),       # gene-transcript1
    re.compile(r'[-_][Pp]\d+$'),                     # gene-P1 (Augustus)
]


def to_gene_id(seq_id: str) -> str:
    """Strip the isoform suffix from a sequence identifier."""
    gene_id = seq_id
    for pat in ISOFORM_PATTERNS:
        stripped = pat.sub("", gene_id)
        if stripped != gene_id:
            return stripped  # apply the first matching pattern and stop
    return gene_id


def keep_longest(input_fasta: str, output_fasta: str) -> tuple[int, int]:
    """
    Returns (n_input, n_output).
    """
    genes: dict = {}

    for rec in SeqIO.parse(input_fasta, "fasta"):
        gene_id = to_gene_id(rec.id)
        if gene_id not in genes or len(rec.seq) > len(genes[gene_id].seq):
            # Re-label the record with the gene ID so downstream tools see
            # a consistent identifier (original full ID preserved in description)
            rec.description = f"original_id={rec.id} {rec.description}"
            rec.id = gene_id
            rec.name = gene_id
            genes[gene_id] = rec

    with open(output_fasta, "w") as out:
        SeqIO.write(genes.values(), out, "fasta")

    return sum(1 for _ in SeqIO.parse(input_fasta, "fasta")), len(genes)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: keep_longest_isoform.py <input.fa> <output.fa>")
        sys.exit(1)

    n_in, n_out = keep_longest(sys.argv[1], sys.argv[2])
    removed = n_in - n_out
    pct = 100 * n_out / max(n_in, 1)
    print(
        f"  {sys.argv[1]:<55}  "
        f"{n_in:>7} isoforms  →  {n_out:>7} genes  "
        f"({removed} isoforms removed, {pct:.1f}% retained)"
    )
PYEOF

# -----------------------------------------------------------------------------
# filter_cestoda_markers.py
# -----------------------------------------------------------------------------
# Selects BUSCO-quality marker genes from the OrthoFinder output.
#
# Selection criteria (adapted from Waterhouse et al. 2007):
#
#   1. PRESENCE: the gene must be detectable in at least MIN_CESTODA_PRESENCE
#      fraction of the 11 Cestoda training species.
#      We relax this from the canonical 90% to 80% because our training
#      assemblies have genuine incompleteness (67–77% BUSCO completeness with
#      platyhelminthes_odb10), so some absences reflect assembly gaps rather
#      than true gene loss.
#
#   2. SINGLE-COPY: among Cestoda species where the gene is present, at least
#      MIN_SINGLECOPY_FRAC must carry exactly one copy.
#      This ensures that the marker is not a member of a tandem duplicate
#      array or a rapidly evolving gene family, which would compromise the
#      one-to-one relationship that BUSCO scoring requires.
#
# Conservation annotation (non-filtering metadata):
#   Each retained marker is tagged based on outgroup presence:
#     "platyhelminthes_conserved" — present in ≥OUTGROUP_CONSERVATION_MIN
#                                   outgroup species; likely ancient gene
#     "cestoda_specific"          — rare or absent in outgroups; potentially
#                                   Cestoda-lineage innovation or highly
#                                   diverged in outgroups
#
# This annotation is written to marker_metadata.tsv and is used during
# interpretation of BUSCO results: low scores on "cestoda_specific" markers
# may reflect genuine Cestoda gene evolution, whereas low scores on
# "platyhelminthes_conserved" markers more likely indicate assembly quality
# problems.
# -----------------------------------------------------------------------------
cat > "${SCRIPTS_DIR}/filter_cestoda_markers.py" << 'PYEOF'
#!/usr/bin/env python3
"""
filter_cestoda_markers.py
Clade-aware BUSCO marker selection from OrthoFinder orthogroups.

Usage:
    python filter_cestoda_markers.py \
        --og_table          <Orthogroups.tsv> \
        --og_seqdir         <Orthogroup_Sequences/> \
        --out_dir           <output_dir> \
        --min_presence      <float> \
        --min_singlecopy    <float> \
        --outgroup_min      <int>
"""

import argparse
import math
import os
import shutil
import sys
from pathlib import Path

import pandas as pd


# ── Clade membership ──────────────────────────────────────────────────────────
# Tags must match the filenames used in 01_retrieve_data.sh (without .fa).
# OrthoFinder uses the input filename (minus extension) as the species label.

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
    """Count gene copies from an OrthoFinder table cell (comma-separated IDs)."""
    if pd.isna(cell):
        return 0
    value = str(cell).strip()
    if value == "" or value.lower() in {"nan", "none"}:
        return 0
    return len([x for x in value.split(",") if x.strip()])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--og_table",        required=True)
    p.add_argument("--og_seqdir",       required=True)
    p.add_argument("--out_dir",         required=True)
    p.add_argument("--min_presence",    type=float, default=0.80)
    p.add_argument("--min_singlecopy",  type=float, default=0.90)
    p.add_argument("--outgroup_min",    type=int,   default=3)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"\n  Reading orthogroup table: {args.og_table}")
    df = pd.read_csv(args.og_table, sep="\t", index_col=0)
    n_total = len(df)
    print(f"  Total orthogroups: {n_total:,}")

    # Map table columns to clade sets
    cestoda_cols   = [c for c in df.columns if c in CESTODA]
    outgroup_cols  = [c for c in df.columns if c in OUTGROUPS]
    unknown_cols   = [c for c in df.columns if c not in CESTODA and c not in OUTGROUPS]

    print(f"\n  Species detected in table:")
    print(f"    Cestoda   : {len(cestoda_cols):2d} / {len(CESTODA)} expected")
    print(f"    Outgroups : {len(outgroup_cols):2d}")

    if unknown_cols:
        print(f"\n  WARNING — unrecognised column(s): {unknown_cols}")
        print("    These will be ignored. Check that column names match species tags.")

    if len(cestoda_cols) == 0:
        print("\n  ERROR — no Cestoda species found in the table.")
        print("    Verify that OrthoFinder input filenames match the CESTODA set.")
        sys.exit(1)

    n_cestoda = len(cestoda_cols)
    presence_min = math.ceil(n_cestoda * args.min_presence)  # absolute count, rounded up

    print(f"\n  Filtering thresholds:")
    print(f"    Min Cestoda presence  : {args.min_presence*100:.0f}% = ≥{presence_min} / {n_cestoda} species")
    print(f"    Min single-copy frac  : {args.min_singlecopy*100:.0f}% of present species")
    print(f"    Outgroup conservation : ≥{args.outgroup_min} outgroup species")

    kept = []
    metadata = []
    removed_presence   = 0
    removed_multicopy  = 0

    for og_id, row in df.iterrows():

        # ── Per-species copy counts (Cestoda only) ────────────────────────────
        copies = {sp: count_copies(row[sp]) for sp in cestoda_cols}
        n_present     = sum(1 for c in copies.values() if c >= 1)
        n_single_copy = sum(1 for c in copies.values() if c == 1)

        # ── Filter 1: presence ────────────────────────────────────────────────
        if n_present < presence_min:
            removed_presence += 1
            continue

        # ── Filter 2: single-copy ─────────────────────────────────────────────
        sc_frac = n_single_copy / n_present
        if sc_frac < args.min_singlecopy:
            removed_multicopy += 1
            continue

        # ── Conservation annotation (does not affect filtering) ───────────────
        n_outgroup = sum(1 for sp in outgroup_cols if count_copies(row[sp]) >= 1)
        conservation = (
            "platyhelminthes_conserved"
            if n_outgroup >= args.outgroup_min
            else "cestoda_specific"
        )

        kept.append(og_id)
        metadata.append({
            "orthogroup_id":        og_id,
            "n_cestoda_present":    n_present,
            "n_cestoda_singlecopy": n_single_copy,
            "singlecopy_fraction":  round(sc_frac, 4),
            "n_outgroup_present":   n_outgroup,
            "conservation":         conservation,
        })

    # ── Summary ───────────────────────────────────────────────────────────────
    n_kept = len(kept)
    n_conserved = sum(1 for m in metadata if m["conservation"] == "platyhelminthes_conserved")
    n_specific  = sum(1 for m in metadata if m["conservation"] == "cestoda_specific")

    print(f"\n  Results:")
    print(f"    Total orthogroups        : {n_total:,}")
    print(f"    Removed (low presence)   : {removed_presence:,}")
    print(f"    Removed (multi-copy)     : {removed_multicopy:,}")
    print(f"    Retained as markers      : {n_kept:,}")
    print(f"")
    print(f"    Platyhelminthes-conserved: {n_conserved:,}")
    print(f"    Cestoda-specific         : {n_specific:,}")

    if n_kept < 100:
        print(f"\n  WARNING: only {n_kept} markers retained.")
        print("    Consider relaxing --min_presence or --min_singlecopy.")
        print("    A functional BUSCO dataset typically needs ≥200 markers.")

    # ── Write metadata ────────────────────────────────────────────────────────
    meta_path = os.path.join(args.out_dir, "marker_metadata.tsv")
    pd.DataFrame(metadata).to_csv(meta_path, sep="\t", index=False)
    print(f"\n  Metadata: {meta_path}")

    # ── Copy FASTA files ──────────────────────────────────────────────────────
    seqdir  = Path(args.og_seqdir)
    copied  = 0
    missing = 0
    for og in kept:
        src = seqdir / f"{og}.fa"
        if src.exists():
            shutil.copy(src, os.path.join(args.out_dir, f"{og}.fa"))
            copied += 1
        else:
            print(f"  WARN: FASTA not found for {og} in {seqdir}")
            missing += 1

    print(f"  FASTAs: {copied} copied, {missing} missing")
    print("  Done.")


if __name__ == "__main__":
    main()
PYEOF

chmod +x "${SCRIPTS_DIR}/clean_sequences.py"
chmod +x "${SCRIPTS_DIR}/keep_longest_isoform.py"
chmod +x "${SCRIPTS_DIR}/filter_cestoda_markers.py"

echo "  Scripts written to ${SCRIPTS_DIR}/"

# =============================================================================
# STAGE 1 — Sequence quality control and cleaning
# =============================================================================
# Remove short sequences, internal stops, high-X sequences, and duplicates.
# Input:  data/01_raw_proteomes/
# Output: data/02_clean_proteomes/
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 1 — Sequence quality control"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

QC_LOG="${LOG_DIR}/qc_summary.tsv"
echo -e "species\tn_raw\tn_clean\tpct_kept\tn_short\tn_internal_stop\tn_high_x\tn_duplicate" > "${QC_LOG}"

for RAW_FA in "${RAW_PROTEOMES}"/*.fa; do
    TAG=$(basename "${RAW_FA}" .fa)
    CLEAN_FA="${CLEAN_PROTEOMES}/${TAG}.fa"
    REPORT="${CLEAN_PROTEOMES}/${TAG}.qc.tsv"

    if [[ -f "${CLEAN_FA}" ]] && [[ -s "${CLEAN_FA}" ]]; then
        echo "  [SKIP] ${TAG}"
        continue
    fi

    python "${SCRIPTS_DIR}/clean_sequences.py" \
        --input      "${RAW_FA}" \
        --output     "${CLEAN_FA}" \
        --report     "${REPORT}" \
        --min_length "${MIN_PROTEIN_LENGTH}" \
        --max_x_frac "${MAX_X_FRACTION}"
done

for REPORT in "${CLEAN_PROTEOMES}"/*.qc.tsv; do
    [[ -f "${REPORT}" ]] || continue
    TAG=$(basename "${REPORT}" .qc.tsv)

    N_RAW=$(awk 'NR>1{n++} END{print n+0}' "${REPORT}")
    N_CLEAN_TAG=$(awk -F'\t' 'NR>1 && $3=="kept"{n++} END{print n+0}' "${REPORT}")
    N_SHORT=$(awk -F'\t' 'NR>1 && $3 ~ /^too_short/{n++} END{print n+0}' "${REPORT}")
    N_STOP=$(awk -F'\t' 'NR>1 && $3 ~ /^internal_stop/{n++} END{print n+0}' "${REPORT}")
    N_X=$(awk -F'\t' 'NR>1 && $3 ~ /^high_X/{n++} END{print n+0}' "${REPORT}")
    N_DUP=$(awk -F'\t' 'NR>1 && $3 == "duplicate"{n++} END{print n+0}' "${REPORT}")
    PCT=$(awk -v kept="${N_CLEAN_TAG}" -v raw="${N_RAW}" 'BEGIN{ if (raw==0) print "0.00"; else printf "%.2f", (100*kept/raw) }')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"         "${TAG}" "${N_RAW}" "${N_CLEAN_TAG}" "${PCT}" "${N_SHORT}" "${N_STOP}" "${N_X}" "${N_DUP}"         >> "${QC_LOG}"
done

N_CLEAN=$(find "${CLEAN_PROTEOMES}" -name "*.fa" -size +0c | wc -l)
echo ""
echo "  Cleaned proteomes: ${N_CLEAN} files in ${CLEAN_PROTEOMES}/"
echo "  QC summary       : ${QC_LOG}"

# =============================================================================
# STAGE 2 — Longest-isoform filtering
# =============================================================================
# Reduce each proteome to one sequence per gene locus.
# Input:  data/02_clean_proteomes/
# Output: data/03_filtered_proteomes/
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 2 — Longest-isoform filtering"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

for CLEAN_FA in "${CLEAN_PROTEOMES}"/*.fa; do
    TAG=$(basename "${CLEAN_FA}" .fa)
    FILT_FA="${FILTERED_PROTEOMES}/${TAG}.fa"

    if [[ -f "${FILT_FA}" ]] && [[ -s "${FILT_FA}" ]]; then
        echo "  [SKIP] ${TAG}"
        continue
    fi

    python "${SCRIPTS_DIR}/keep_longest_isoform.py" "${CLEAN_FA}" "${FILT_FA}"
done

N_FILT=$(find "${FILTERED_PROTEOMES}" -name "*.fa" -size +0c | wc -l)
echo ""
echo "  Filtered proteomes: ${N_FILT} files in ${FILTERED_PROTEOMES}/"

# Final sequence count per species
echo ""
echo "  Final sequence counts (after QC + isoform filtering):"
printf "  %-50s  %s\n" "Species" "Sequences"
printf "  %-50s  %s\n" "-------" "---------"
for FA in "${FILTERED_PROTEOMES}"/*.fa; do
    [[ -f "${FA}" ]] || continue
    printf "  %-50s  %d\n" "$(basename "${FA}" .fa)" "$(grep -c '^>' "${FA}")"
done

# =============================================================================
# STAGE 3 — OrthoFinder
# =============================================================================
# OrthoFinder performs:
#   1. All-vs-all DIAMOND protein search (ultra-sensitive mode)
#   2. Normalisation of bit-scores to correct for protein length / species bias
#   3. Orthogroup inference via the MCL graph-clustering algorithm
#   4. Multiple sequence alignment per orthogroup (MAFFT; -M msa flag)
#   5. Gene trees and rooted species tree (FastTree)
#
# Key parameters:
#   -S diamond_ultra_sens  Maximises sensitivity for divergent flatworm proteins
#                          at modest speed cost (still faster than BLAST).
#   -M msa                 MSA-based tree inference (more accurate than default
#                          UPGMA-based; required for rooted species tree).
#   -A mafft               MAFFT for per-orthogroup alignments within OrthoFinder.
#   -T fasttree            FastTree for gene-tree inference.
#   -t                     Threads for DIAMOND and MAFFT parallelisation.
#   -a                     Threads for parallel OrthoFinder analysis steps.
#   -n                     Run label (determines output subdirectory name).
#
# OrthoFinder is resumable: if a previous run exists with the same -n label,
# delete the Results_* directory and rerun to restart from scratch.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 3 — OrthoFinder"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

OF_RUN_NAME="cestoda_run1"
OF_RESULTS_DIR="${OG_DIR}/Results_${OF_RUN_NAME}"

if [[ -d "${OF_RESULTS_DIR}" ]]; then
    echo "  [SKIP] OrthoFinder results already exist: ${OF_RESULTS_DIR}"
    echo "         Delete that directory and rerun to restart OrthoFinder."
else
    echo "  Starting OrthoFinder..."
    echo "  Threads: ${OF_THREADS} (DIAMOND/MAFFT) × ${OF_ANALYSIS_THREADS} (analysis)"
    echo "  Input:   ${FILTERED_PROTEOMES}/"
    echo "  Estimated runtime: 6–20 hours (depends on CPU count and proteome sizes)"
    echo ""

    # OrthoFinder refuses a user-specified -o directory if it already exists,
    # even when it is empty. Because 00_setup_environment.sh creates the whole
    # project tree in advance, we must remove an empty or stale OG_DIR before
    # launching a fresh run. We only auto-remove it when there is no completed
    # Results_* directory inside; otherwise we leave the completed run intact.
    if [[ -d "${OG_DIR}" ]]; then
        if find "${OG_DIR}" -mindepth 1 -maxdepth 1 -type d -name 'Results_*' | grep -q .; then
            echo "  [ERROR] ${OG_DIR} already contains Results_* directories, but ${OF_RESULTS_DIR} was not found."
            echo "          Inspect ${OG_DIR} and either:"
            echo "            1) point the script to the correct results directory, or"
            echo "            2) remove ${OG_DIR} to restart OrthoFinder from scratch."
            exit 1
        fi

        echo "  Removing pre-existing OrthoFinder output directory so -o can be created fresh:"
        echo "    ${OG_DIR}"
        rm -rf "${OG_DIR}"
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
        -o "${OG_DIR}" \
        2>&1 | tee "${LOG_DIR}/orthofinder.log"

    echo ""

    if [[ ! -d "${OF_RESULTS_DIR}" ]]; then
        echo "  [ERROR] OrthoFinder did not produce the expected results directory:"
        echo "          ${OF_RESULTS_DIR}"
        echo "          Check ${LOG_DIR}/orthofinder.log"
        exit 1
    fi

    echo "  OrthoFinder complete."
fi

# Locate the key output files — OrthoFinder may place them in subdirectories
# whose names depend on the version and run parameters.
OG_TABLE=$(find "${OF_RESULTS_DIR}" -name "Orthogroups.tsv"         -not -path "*/__pycache__/*" | head -1)
OG_SEQDIR=$(find "${OF_RESULTS_DIR}" -name "Orthogroup_Sequences"   -type d | head -1)

if [[ -z "${OG_TABLE}" ]]; then
    echo ""
    echo "  [ERROR] Orthogroups.tsv not found under ${OF_RESULTS_DIR}/"
    echo "          Check OrthoFinder log: ${LOG_DIR}/orthofinder.log"
    exit 1
fi

if [[ -z "${OG_SEQDIR}" ]]; then
    echo ""
    echo "  [ERROR] Orthogroup_Sequences/ not found under ${OF_RESULTS_DIR}/"
    echo "          OrthoFinder may not have completed successfully."
    exit 1
fi

N_ORTHOGROUPS=$(($(wc -l < "${OG_TABLE}") - 1))  # subtract header line
echo ""
echo "  Orthogroup table : ${OG_TABLE}"
echo "  Sequence dir     : ${OG_SEQDIR}"
printf "  Total orthogroups: %s\n" "${N_ORTHOGROUPS}"

# =============================================================================
# STAGE 4 — Marker selection
# =============================================================================
# Filter orthogroups to retain BUSCO-quality Cestoda marker genes.
# Input:  Orthogroups.tsv + Orthogroup_Sequences/
# Output: data/05_cestoda_markers/*.fa + marker_metadata.tsv
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 4 — Marker selection"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

MARKER_META="${MARKER_DIR}/marker_metadata.tsv"
EXISTING_MARKERS=$(find "${MARKER_DIR}" -name "OG*.fa" -size +0c | wc -l)

if [[ -s "${MARKER_META}" ]] && [[ "${EXISTING_MARKERS}" -gt 0 ]]; then
    echo "  [SKIP] Marker FASTAs already exist in ${MARKER_DIR}/"
    echo "         Metadata: ${MARKER_META}"
else
    # Clean partial outputs from a previous interrupted run before regenerating.
    rm -f "${MARKER_DIR}"/OG*.fa "${MARKER_DIR}/marker_metadata.tsv"

    python "${SCRIPTS_DIR}/filter_cestoda_markers.py" \
        --og_table       "${OG_TABLE}" \
        --og_seqdir      "${OG_SEQDIR}" \
        --out_dir        "${MARKER_DIR}" \
        --min_presence   "${MIN_CESTODA_PRESENCE}" \
        --min_singlecopy "${MIN_SINGLECOPY_FRAC}" \
        --outgroup_min   "${OUTGROUP_CONSERVATION_MIN}" \
        2>&1 | tee "${LOG_DIR}/filter_markers.log"
fi

N_MARKERS=$(find "${MARKER_DIR}" -name "OG*.fa" -size +0c | wc -l)
echo ""
echo "  Markers retained: ${N_MARKERS}"

if [[ "${N_MARKERS}" -eq 0 ]]; then
    echo "  [ERROR] No marker FASTAs produced."
    echo "          Check ${LOG_DIR}/filter_markers.log"
    exit 1
fi

# =============================================================================
# STAGE 5 — Multiple sequence alignment (MAFFT)
# =============================================================================
# Re-align each marker FASTA using MAFFT --auto --reorder.
#
# We use a bounded parallel pool instead of `--thread -1` so the total thread
# count stays under explicit control. Each MAFFT job gets THREADS_PER_JOB threads
# and we run at most PARALLEL_JOBS jobs concurrently.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 5 — Multiple sequence alignment (MAFFT)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Jobs: ${PARALLEL_JOBS} parallel × ${THREADS_PER_JOB} threads each"
echo ""

run_mafft_pool() {
    local PIDS=()
    local N_ALIGN_NEW=0
    local N_ALIGN_SKIP=0
    local N_ALIGN_FAIL=0

    for MARKER_FA in "${MARKER_DIR}"/OG*.fa; do
        local OG
        local MSA
        local N_SEQS

        OG=$(basename "${MARKER_FA}" .fa)
        MSA="${ALIGN_DIR}/${OG}.msa"

        if [[ -f "${MSA}" ]] && [[ -s "${MSA}" ]]; then
            N_ALIGN_SKIP=$((N_ALIGN_SKIP + 1))
            continue
        fi

        N_SEQS=$(grep -c '^>' "${MARKER_FA}" 2>/dev/null || echo 0)
        if [[ "${N_SEQS}" -lt 2 ]]; then
            echo "  [SKIP] ${OG}: only ${N_SEQS} sequence(s) — need ≥2 for alignment"
            N_ALIGN_FAIL=$((N_ALIGN_FAIL + 1))
            continue
        fi

        (
            mafft                 --auto                 --reorder                 --quiet                 --thread "${THREADS_PER_JOB}"                 "${MARKER_FA}"                 > "${MSA}" 2>>"${LOG_DIR}/mafft.log"

            [[ -s "${MSA}" ]]
        ) &
        PIDS+=($!)
        N_ALIGN_NEW=$((N_ALIGN_NEW + 1))

        if [[ "${#PIDS[@]}" -ge "${PARALLEL_JOBS}" ]]; then
            for PID in "${PIDS[@]}"; do
                if ! wait "${PID}"; then
                    N_ALIGN_FAIL=$((N_ALIGN_FAIL + 1))
                fi
            done
            PIDS=()
            echo "  Aligned ${N_ALIGN_NEW} marker families so far..."
        fi
    done

    for PID in "${PIDS[@]}"; do
        if ! wait "${PID}"; then
            N_ALIGN_FAIL=$((N_ALIGN_FAIL + 1))
        fi
    done

    N_ALIGN_TOTAL=$(find "${ALIGN_DIR}" -name "*.msa" -size +0c | wc -l)
    echo "  Alignments: ${N_ALIGN_NEW} new / ${N_ALIGN_SKIP} existing / ${N_ALIGN_FAIL} failed-or-skipped"
    echo "  Total MSA files: ${N_ALIGN_TOTAL}"
}

run_mafft_pool

# =============================================================================
# STAGE 6 — HMM profile construction (HMMER hmmbuild)
# =============================================================================
# hmmbuild reads a multiple sequence alignment and constructs a profile HMM
# that captures the position-specific amino acid frequencies and gap patterns
# of the marker gene family.
#
# These profiles are the core of BUSCO scoring: when BUSCO searches a new
# genome, it uses hmmsearch with these profiles to identify candidate gene
# copies and assess their quality.
#
# Options:
#   --cpu N    : threads per hmmbuild instance
#   --amino    : explicit protein mode (avoids auto-detection edge cases)
#   -n <name>  : label the profile with the orthogroup ID
#
# Parallelisation:
#   We run PARALLEL_JOBS instances of hmmbuild concurrently, each using
#   THREADS_PER_JOB threads, for a total of PARALLEL_JOBS × THREADS_PER_JOB
#   active threads.  Adjust these values in project.cfg to fit your system.
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stage 6 — HMM profile construction (hmmbuild)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Jobs: ${PARALLEL_JOBS} parallel × ${THREADS_PER_JOB} threads each"
echo ""

# Simple parallel job pool: maintain at most PARALLEL_JOBS background PIDs.
run_parallel_jobs() {
    local PIDS=()
    local N_BUILT=0
    local N_SKIP=0

    for MSA in "${ALIGN_DIR}"/*.msa; do
        OG=$(basename "${MSA}" .msa)
        HMM="${HMM_DIR}/${OG}.hmm"

        if [[ -f "${HMM}" ]] && [[ -s "${HMM}" ]]; then
            N_SKIP=$((N_SKIP + 1))
            continue
        fi

        # Spawn a background hmmbuild job
        hmmbuild \
            --cpu "${THREADS_PER_JOB}" \
            --amino \
            -n "${OG}" \
            "${HMM}" \
            "${MSA}" \
            >> "${LOG_DIR}/hmmbuild.log" 2>&1 &

        PIDS+=($!)
        N_BUILT=$((N_BUILT + 1))

        # When the pool is full, wait for all current jobs before launching more
        if [[ "${#PIDS[@]}" -ge "${PARALLEL_JOBS}" ]]; then
            for PID in "${PIDS[@]}"; do
                wait "${PID}" || {
                    echo "  [WARN] hmmbuild PID ${PID} exited non-zero — check ${LOG_DIR}/hmmbuild.log"
                }
            done
            PIDS=()
            echo "  Built ${N_BUILT} profiles so far..."
        fi
    done

    # Drain remaining background jobs
    for PID in "${PIDS[@]}"; do
        wait "${PID}" || {
            echo "  [WARN] hmmbuild PID ${PID} exited non-zero"
        }
    done

    echo "  hmmbuild: ${N_BUILT} new profiles built / ${N_SKIP} already existed"
}

run_parallel_jobs

N_HMM_TOTAL=$(find "${HMM_DIR}" -name "*.hmm" -size +0c | wc -l)
echo "  Total HMM profiles: ${N_HMM_TOTAL}"

# =============================================================================
# STAGE SUMMARY
# =============================================================================

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Analysis pipeline complete — $(date)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Species processed     : ${N_FILT}"
echo "  Orthogroups inferred  : ${N_ORTHOGROUPS}"
echo "  Markers selected      : ${N_MARKERS}"
echo "  Alignments produced   : ${N_ALIGN_TOTAL}"
echo "  HMM profiles built    : ${N_HMM_TOTAL}"
echo ""
echo "  Log: ${LOG}"
echo "  Next step: bash 03_build_odb.sh"
echo ""
