# Cestoda BUSCO Custom Lineage Builder — Full Pipeline

## Overview

This repository contains a **complete, production-grade pipeline** to build a **custom BUSCO v6 lineage dataset for Cestoda**, from raw proteomes to validated BUSCO-compatible lineage.

The pipeline is designed to ensure:

- **Strict identifier contract consistency**
- **Reproducibility and resumability**
- **Robust error handling**
- **Functional validity in BUSCO (not just structural correctness)**

---

## What This Pipeline Solves

Standard custom BUSCO lineage attempts often fail with:

- “No jobs to run on hmmsearch”
- 0 complete BUSCOs
- 100% missing BUSCOs

These failures are caused by **identifier contract inconsistencies**, including mismatches between:

- FASTA IDs
- Orthogroup IDs
- HMM filenames
- HMM internal NAME fields
- scores_cutoff
- refseq_db
- BUSCO internal expectations

This pipeline enforces a **single consistent identifier system end-to-end**.

---

## Full Pipeline Structure

```
00_setup_environment.sh
01_retrieve_data.sh
02_run_analysis.sh
03_build_odb.sh
04_test_cestoda_lineage.sh
```

### Step Descriptions

#### 00_setup_environment.sh
- Creates conda environment
- Installs required tools:
  - BUSCO v6
  - HMMER
  - OrthoFinder
  - Python dependencies

#### 01_retrieve_data.sh
- Downloads proteomes used for lineage construction
- Organizes input data

#### 02_run_analysis.sh
Core analytical pipeline:

- Sequence QC
- Longest isoform filtering
- FASTA identifier audit (critical checkpoint)
- OrthoFinder clustering
- Alignment (MAFFT)
- HMM construction
- hmmsearch execution

#### 03_build_odb.sh
Lineage export:

- Enforces **BUSCO identifier contract**
- Generates:
  - HMMs (renamed + NAME fixed)
  - scores_cutoff
  - refseq_db.faa.gz
  - ancestral sequences
  - metadata files
- **Outgroup-aware filtering applied**
- Full contract validation step

#### 04_test_cestoda_lineage_both_modes.sh
Validation:

- Runs BUSCO in:
  - protein mode
  - genome mode
- Validates:
  - dataset functionality
  - performance metrics
  - parsing robustness

---

## Key Design Principles

### 1. Identifier Contract Integrity

A single BUSCO ID is propagated across:

- HMM filenames
- HMM NAME fields
- scores_cutoff
- ogs.id.info
- links_to_ODB12.txt
- refseq_db headers

No mixing of OG IDs and BUSCO IDs is allowed.

---

### 2. Outgroup-Aware Marker Filtering

Outgroups are used explicitly to improve lineage specificity.

Filtering criteria:

```
max_outgroup_presence_frac = 0.40
max_outgroup_singlecopy_frac = 0.50
```

Effect:

- Removes overly conserved genes
- Retains lineage-informative markers

---

### 3. Strict FASTA Validation

Pipeline fails early if:

- Duplicate IDs exist
- Empty IDs exist
- Illegal characters detected

---

### 4. No Silent Reuse of Broken Outputs

Critical steps reset outputs when needed to avoid:

- stale artifacts
- hidden inconsistencies

---

### 5. Functional Validation (Not Just Structure)

Final dataset must:

- Run in BUSCO genome mode
- Run in BUSCO protein mode
- Produce biologically meaningful scores

---

## How to Run

### 1. Setup environment

```bash
bash 00_setup_environment.sh
```

### 2. Retrieve data

```bash
bash 01_retrieve_data.sh
```

### 3. Run analysis

```bash
bash 02_run_analysis.sh
```

### 4. Build lineage

```bash
bash 03_build_odb.sh
```

### 5. Validate lineage

```bash
bash 04_test_cestoda_lineage.sh
```

---

## Output

Final dataset:

```
cestoda_odb_custom_odb12/
```

Contains:

- hmms/
- scores_cutoff
- refseq_db.faa.gz
- ancestral
- dataset.cfg
- info/

---

## Validation Results (Example)

### Genome mode

- ~54–85% Complete BUSCOs (depending on species)
- Functional hmmsearch execution

### Protein mode

- ~80–85% Complete BUSCOs
- Low missing fraction

---

## Important Notes

- Step 03 can be rerun without repeating Step 02
- Heavy computations are only in Step 02
- Pipeline is designed for HPC environments

---

## Requirements

- Linux environment
- Conda
- BUSCO v6
- HMMER
- OrthoFinder
- MAFFT
- Python 3

---

## License

MIT License

---

## Author

Dr. Lucas L. Maldonado
