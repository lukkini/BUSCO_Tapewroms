# Cestoda BUSCO Custom Lineage Builder (with Outgroup Filtering)

## Overview

This repository provides a **production-grade pipeline** to build a custom BUSCO lineage dataset for **Cestoda**, including:

- Orthogroup-based marker selection
- HMM construction and validation
- **Explicit outgroup-aware filtering**
- BUSCO v6-compatible dataset export
- Full contract validation

This version ensures that **outgroups are treated as true outgroups**, not just additional taxa.

---

## Key Feature: Outgroup-Aware Filtering

Before exporting the final lineage (`cestoda_odb_custom_odb12`), markers are filtered using:

- Outgroup presence constraint
- Outgroup single-copy behavior

### Filtering criteria

```text
max_outgroup_presence_frac = 0.40
max_outgroup_singlecopy_frac = 0.50
```

This ensures:
- Markers are conserved in Cestoda
- Markers are not broadly conserved across outgroups
- Lineage specificity is enforced

---

## Pipeline Structure

```
01_retrieve_data.sh
02_run_analysis.sh
03_build_odb_outgroup_filtered.sh   <-- THIS SCRIPT
04_test_cestoda_lineage.sh
```

---

## Download Script

Download the build script directly:

👉 [Download 03_build_odb_outgroup_filtered.sh](./03_build_odb_outgroup_filtered.sh)

Or via command line:

```bash
wget https://raw.githubusercontent.com/<YOUR_REPO>/main/03_build_odb_outgroup_filtered.sh
chmod +x 03_build_odb_outgroup_filtered.sh
```

---

## Usage

Run the lineage build:

```bash
bash 03_build_odb_outgroup_filtered.sh
```

---

## Requirements

- BUSCO v6
- HMMER
- Python 3
- OrthoFinder (already used in previous steps)
- hmmemit (from HMMER)

---

## Output

The final dataset will be created at:

```
cestoda_odb_custom_odb12/
```

Containing:

- `hmms/` → HMM profiles (BUSCO ID renamed)
- `scores_cutoff`
- `links_to_ODB12.txt`
- `refseq_db.faa.gz`
- `ancestral`
- `dataset.cfg`
- `info/`
  - `busco_id_map.tsv`
  - `outgroup_filter_report.tsv`
  - `ogs.id.info`

---

## Validation

After building, run:

```bash
bash 04_test_cestoda_lineage.sh
```

---

## Notes

- This script **does NOT recompute upstream steps**
- It only rebuilds the **final lineage export**
- Safe to run without re-running step 02

---

## Rationale

Standard pipelines often:
- include outgroups during orthogroup construction
- but **do not enforce outgroup exclusion at export time**

This implementation fixes that by:
- explicitly filtering markers using outgroup signal
- ensuring lineage specificity

---

## License

MIT License

---

## Author

Dr. Lucas L. Maldonado
