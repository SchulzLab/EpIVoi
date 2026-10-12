# Chromosome 21 example dataset

This directory contains the input files for the chromosome 21
EpIVoi example based on IHEC CRE-RF models.

## Files

The original object is:

`epivoi_obj_cre_rf_IHEC_1MB_test_performance_chr21.rds`

It contains model metadata, gene annotations, model performance,
and file indices for 346 genes.

Additional input files are provided in:

- `feature_importance/`: 346 gene-specific feature importance files
- `epigenetic_signal/`: 346 gene-specific epigenetic signal files

## Prepare the example

From the repository root, run:

```bash
Rscript scripts/prepare_chr21_example.R
This creates epivoi_object/data/epivoi_chr21_local.rds.
The script updates the file references to relative paths without
modifying the original RDS.
The locally generated RDS does not need to be committed to Git.
Configure external resources
Copy .Renviron.example to .Renviron and configure the
required external resources.
IGV.js assets, RefSeq annotations, the matching hg38 reference
genome and TF analysis dependencies must be available locally.
Start EpIVoi
From the repository root:
Rscript app.R epivoi_object/data/epivoi_chr21_local.rds


The application loads feature importance and epigenetic signal
files from this directory.
TRAP, PASTAA and FIMO require their respective executables,
motif resources and reference genome to be configured separately.
