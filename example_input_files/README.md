# EpIVoi example input dataset

This directory provides a small example dataset for starting and
exploring EpIVoi without downloading the original large model outputs.

The dataset contains nine gene models and two biological groups:
Astrocytes:Alzheimer's and Astrocytes:Unaffected.

## Included files

| File | Description |
| --- | --- |
| `example_epivoi_object.rds` | EpIVoi input object containing model performance, gene annotations and file indices |
| `feature_importance/` | Nine gene-specific feature importance matrices |
| `epigenetic_signal/` | Nine gene-specific epigenetic signal matrices |
| `example_ncbiRefSeq.txt.gz` | RefSeq annotation for the example genes |
| `example_genome_sparse.fa.gz` | Coordinate-preserving, restricted hg38 reference |
| `example_regions.bed` | Genomic regions included in the example |
| `example_energy_matrix.txt` | RUNX1 TRAP binding energy matrix |
| `example_motifs.meme` | RUNX1 motif for FIMO |
| `example_fimo_background.txt` | Example nucleotide background frequencies |
| `example_regulatory_build.bed` | Illustrative regulatory annotation |

The reference FASTA contains the original hg38 sequences at the
selected genomic coordinates. Unused positions are filled with N.
It is not a replacement for the complete hg38 reference genome.

The example is intended for testing and exploring EpIVoi.
Enrichment results are not intended as biological findings.

## Requirements

Install the R dependencies defined by the repository environment,
bedtools, samtools and MEME Suite.

The TRAP and PASTAA executables must also be installed.

IGV.js JavaScript and CSS assets must be available in a local
directory. Their location is configured through EPIVOI_IGV_ASSETS_DIR.

## Prepare the reference

From the repository root, run:

```bash
gzip -dkf example_input_files/example_genome_sparse.fa.gz
samtools faidx example_input_files/example_genome_sparse.fa
```

The uncompressed reference occupies approximately 540 MB.
## Configure resources

Set the following environment variables in your local `.Renviron` or shell:

```dotenv
EPIVOI_OBJECT_PATH=/absolute/path/to/EpIVoi/example_input_files/example_epivoi_object.rds
EPIVOI_REFSEQ_FILE=/absolute/path/to/EpIVoi/example_input_files/example_ncbiRefSeq.txt.gz
EPIVOI_GENOME_FASTA=/absolute/path/to/EpIVoi/example_input_files/example_genome_sparse.fa
EPIVOI_PASTAA_ENERGY_MATRIX=/absolute/path/to/EpIVoi/example_input_files/example_energy_matrix.txt
EPIVOI_MEME_MOTIF_FILE=/absolute/path/to/EpIVoi/example_input_files/example_motifs.meme
EPIVOI_FIMO_ENCODE_BG=/absolute/path/to/EpIVoi/example_input_files/example_fimo_background.txt

```

Configure `EPIVOI_TRAP_BIN`, EPIVOI_PASTAA_BIN and
EPIVOI_IGV_ASSETS_DIR for your local installation.
EPIVOI_REGULATORY_BUILD_BED is optional. When no local file is
configured, EpIVoi can retrieve supported Ensembl Regulatory
Build annotations dynamically.
## Start the application

From the repository root, run:

```bash
Rscript app.R example_input_files/example_epivoi_object.rds
```


The example object references its feature importance and epigenetic
signal files using paths relative to the object location.
## Available analyses
The example can be used to explore model performance, feature
importance, genomic regions and transcription factor enrichment.
The small motif collection contains RUNX1 only. PASTAA may
complete successfully without finding statistically significant
transcription factors.
The example data are intended to demonstrate application
functionality rather than support biological interpretation.
