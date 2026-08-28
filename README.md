# EpIVoi

EpIVoi is an interactive **R/Shiny application** for exploring model performance, feature importance, epigenetic signal, transcription-factor enrichment, and genomic context from EpIVoi visualization objects.

The application combines model-level and gene-level summaries with interactive genomic visualization and downstream TF analysis.

![EpIVoi overview](docs/screenshots/01_overview.pdf)

## Main features

- Gene-level model performance overview
- Train-versus-test performance comparison
- Feature Importance / SHAP region visualization
- Interactive IGV genomic view
- Epigenetic-signal tracks
- Ensembl regulatory annotations
- PASTAA/TRAP transcription-factor enrichment
- Gene-level and biological-group-level TF enrichment
- Positive, negative, and absolute Feature Importance analyses
- Combined positive/negative TF-enrichment dot plots
- FIMO motif-site prediction and IGV visualization
- Downloadable plots, tables, and reproducibility code

---

## Application overview

EpIVoi is organized into several tabs. Each tab addresses a different part of the analysis workflow.

### 1. Overview

The **Overview** tab summarizes model performance across all genes.

Users can:

- inspect the distribution of a selected performance metric;
- switch between density and histogram views;
- compare performance across biological groups;
- classify genes into `failed`, `low`, `medium`, and `high` performance classes;
- adjust the correlation cutoffs used for these classes;
- search for a gene by gene symbol or Ensembl ID;
- inspect and download the performance table and plots.

This tab is useful for identifying well-performing models before inspecting individual genes in more detail.

![Overview tab](docs/screenshots/01_overview.pdf)

### 2. Train vs Test

The **Train vs Test** tab compares model performance on the training and test datasets across genes.

Users can:

- select the train/test metric;
- compare train and test performance in a scatter plot;
- optionally use a log scale for error metrics;
- filter points by test-correlation class;
- highlight the top-performing genes;
- inspect and download the corresponding gene table.

This view can help identify differences between training and test performance and select genes for downstream inspection.

![Train vs Test tab](docs/screenshots/02_train_test.pdf)

### 3. Feature Importance

The **Feature Importance** tab visualizes the genomic regions that contribute most strongly to the model prediction for a selected gene and biological group.

Users can:

- select a biological group;
- search for a gene;
- choose the number of top regions to display;
- filter regions by minimum absolute Feature Importance;
- inspect signed Feature Importance values in the plot and table;
- download the selected regions and visualization.

Regions are ranked by the absolute magnitude of Feature Importance while the sign of the original score is retained.

![Feature Importance tab](docs/screenshots/03_feature_importance.pdf)

### 4. IGV

The **IGV** tab places the selected gene and model-derived regions into their genomic context.

The view can include:

- gene annotation;
- SHAP / Feature Importance tracks;
- epigenetic-signal tracks such as ATAC;
- Ensembl regulatory annotations;
- FIMO-predicted TF-binding sites.

Users can select one or more biological groups, load a gene locus, optionally display the Ensembl regulatory track, run FIMO, and export selected genomic tracks as a vector PDF.

A local regulatory-build BED file is not required for the standard configuration. If no local override is supplied, compatible Ensembl regulatory annotations can be obtained dynamically.

![IGV tab](docs/screenshots/04_igv.pdf)

### 5. TF enrichment

The **TF enrichment** tab runs PASTAA/TRAP on genomic regions ranked by Feature Importance.

Two analysis modes are available:

**Run per gene**

A selected gene is analyzed separately for one or more biological groups.

**Run per biological group**

Regions are pooled across genes for the selected biological group. Before pooling, raw Feature Importance values are transformed to within-gene z-scores so that genes with different score scales can be compared more fairly. Gene models can optionally be filtered using a selected test-correlation metric.

The analysis supports:

- positive Feature Importance;
- negative Feature Importance;
- positive and negative analyses;
- absolute Feature Importance.

Positive and negative TF enrichment can be compared in a **combined dot plot**, while the currently selected direction is also displayed as a TF bar plot. Absolute Feature Importance results are shown in a separate dot plot.

For biological-group analyses, gene-level Feature Importance files can be read in parallel on Unix-like systems.

![TF enrichment tab](docs/screenshots/05_tf_enrichment.pdf)

### 6. FIMO motif analysis

FIMO is launched from the **IGV** workflow.

TFs can be obtained from the latest saved PASTAA results or entered manually. EpIVoi maps TF names to available motifs, runs FIMO using the configured motif database and background, and adds successful motif-site predictions to the IGV view.

![FIMO results](docs/screenshots/06_fimo.pdf)

### 7. Reproducibility

The **Reproducibility** tab provides R code corresponding to the current EpIVoi analysis settings.

The generated code documents the major analysis steps outside the GUI, including Feature Importance loading, biological-group z-score calculation, region ranking, sequence preparation, TRAP/PASTAA analysis, and multiple-testing correction.

The reproducibility script can be downloaded directly from the application.

### 8. Help

The **Help** tab contains an in-app explanation of the main EpIVoi concepts, plots, analysis modes, and interpretation of the results.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/SchulzLab/EpIVoi.git
cd EpIVoi
```

### 2. Create the Conda environment

The repository contains an `environment.yml` file with the required software dependencies.

```bash
conda env create -f environment.yml
conda activate epivoi
```

The environment contains R and the R packages required by the Shiny application, together with MEME Suite utilities used by FIMO.

PASTAA and TRAP are external executables and must be available separately.

### 3. Configure external resources

EpIVoi can use project-local resources or paths supplied through environment variables.

| Variable | Purpose |
|---|---|
| `EPIVOI_OBJECT_PATH` | EpIVoi visualization object (`.rds`) |
| `EPIVOI_REFSEQ_FILE` | RefSeq annotation used for vector/PDF export |
| `EPIVOI_IGV_ASSETS_DIR` | Directory containing `igv.min.js` and `igv.min.css` |
| `EPIVOI_TRAP_BIN` | TRAP executable |
| `EPIVOI_PASTAA_BIN` | PASTAA executable |
| `EPIVOI_PASTAA_ENERGY_MATRIX` | PASTAA energy matrix |
| `EPIVOI_GENOME_FASTA` | Reference-genome FASTA |
| `EPIVOI_FIMO_BIN` | FIMO executable; otherwise searched on `PATH` |
| `EPIVOI_FASTA_GET_MARKOV_BIN` | `fasta-get-markov`; otherwise searched on `PATH` |
| `EPIVOI_MEME_MOTIF_FILE` | MEME-format motif database |
| `EPIVOI_FIMO_ENCODE_BG` | ENCODE background file used by FIMO |
| `EPIVOI_PASTAA_CORES` | Optional worker count for biological-group PASTAA preprocessing |
| `EPIVOI_REGULATORY_BUILD_BED` | Optional local regulatory annotation override |

Example configuration:

```bash
export EPIVOI_OBJECT_PATH="/path/to/epivoi_object.rds"
export EPIVOI_REFSEQ_FILE="/path/to/ncbiRefSeq.txt.gz"
export EPIVOI_IGV_ASSETS_DIR="/path/to/igv"

export EPIVOI_TRAP_BIN="/path/to/TRAP"
export EPIVOI_PASTAA_BIN="/path/to/PASTAA"
export EPIVOI_PASTAA_ENERGY_MATRIX="/path/to/Jaspar_Hocomoco_Kellis_human_energy.txt"

export EPIVOI_GENOME_FASTA="/path/to/hg38.fa"
export EPIVOI_MEME_MOTIF_FILE="/path/to/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.meme"
export EPIVOI_FIMO_ENCODE_BG="/path/to/fimo_encode_background.txt"
```

Optionally set the number of parallel workers used during biological-group PASTAA preprocessing:

```bash
export EPIVOI_PASTAA_CORES=4
```

`EPIVOI_REGULATORY_BUILD_BED` does **not** need to be set for the standard configuration.

### 4. Start EpIVoi

```bash
Rscript app.R
```

Open the Shiny URL printed in the terminal.

---

## Toy example

A complete walkthrough is provided in:

[`docs/toy_example.md`](docs/toy_example.md)

A typical workflow is:

1. select a gene such as `ABHD5`;
2. inspect its model performance in **Overview** and **Train vs Test**;
3. inspect its top Feature Importance regions;
4. load the gene and selected biological group in **IGV**;
5. optionally display Ensembl regulatory annotations;
6. run PASTAA for positive and negative Feature Importance;
7. inspect the combined positive/negative TF-enrichment dot plot;
8. run absolute Feature Importance analysis if desired;
9. use saved PASTAA TFs or manual TF names for FIMO;
10. inspect predicted TF-binding sites in IGV;
11. download plots, tables, and reproducibility code.

---

## Input object

EpIVoi expects an RDS visualization object containing the data required by the enabled application modules.

Depending on the object, these data may include:

- species and genome information;
- gene annotation;
- model-performance data;
- Feature Importance data;
- epigenetic-signal data;
- file indexes and metadata used to locate genomic tracks.

The exact set of available biological groups and tracks therefore depends on the loaded object.

---

## Regulatory annotations

A local regulatory-build BED file is optional.

When `EPIVOI_REGULATORY_BUILD_BED` is not supplied, EpIVoi can retrieve compatible Ensembl regulatory annotations dynamically for supported genomes.

A local file can still be supplied when an offline or fixed regulatory annotation source is preferred.

---

## Parallel biological-group analysis

For biological-group TF enrichment, EpIVoi may need to read Feature Importance data for many genes.

On Unix-like systems this step can be parallelized. The worker count can be set explicitly:

```bash
export EPIVOI_PASTAA_CORES=4
```

If no value is provided, the application selects a suitable number of workers automatically. A sequential fallback is used where fork-based parallelization is unavailable.

---

## Reproducibility and downloads

EpIVoi provides downloads for relevant tables and plots throughout the application.

The **Reproducibility** tab additionally generates an R script reflecting the current analysis settings and major processing steps.

Large external resources are referenced through configurable paths rather than hard-coded into the repository.

---

## Documentation screenshots

Screenshots used by this README are stored in:

```text
docs/screenshots/
```

Expected filenames:

```text
01_overview.pdf
02_train_test.pdf
03_feature_importance.pdf
04_igv.pdf
05_tf_enrichment.pdf
06_fimo.pdf
```

---

## License

Add the project license here if/when a license is selected.
