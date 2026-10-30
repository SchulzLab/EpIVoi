# EpIVoi

EpIVoi is an interactive R/Shiny application for interpreting gene specific models that predict gene expression from epigenomic features. It brings together model performance, genomic region feature importance, epigenetic signal, genome browser tracks and transcription factor (TF) enrichment in a standardized, model agnostic visualization framework.

The application is designed for exploration and hypothesis generation. It does not train prediction models, and feature importance should not be interpreted as evidence of causal regulatory activity.

The main modules are **Performance Overview**, **Feature Importance**, **Genome Browser** and **Transcription Factor Enrichment**, as presented in the accompanying EpIVoi manuscript. The interface also includes a dedicated Train vs Test view, Reproducibility and Help.

## Docker Quick Start

The development branch includes a ready-to-run Docker example with nine genes.

```bash
git clone https://github.com/SchulzLab/EpIVoi.git
cd EpIVoi
git checkout aleksandra-epivoi-updates
docker compose up --build
```

Open http://localhost:3838 in your browser.

No manual preparation of Docker inputs is required.
See [DOCKER_SETUP.md](DOCKER_SETUP.md) for detailed instructions.

The following instructions describe installation outside Docker.

## Installation

### Clone the repository

```bash
git clone https://github.com/SchulzLab/EpIVoi.git
cd EpIVoi
```

To use the current development branch containing the chromosome 21 example files and TF enrichment updates, run:

```bash
git checkout aleksandra-epivoi-updates
```

### Install the environment

```bash
conda env create -f environment.yml
conda activate epivoi
```

The `environment.yml` describes the software environment. Install `bedtools` as well, because the TF workflows call `bedtools getfasta`. For local installation outside Docker, TRAP and PASTAA must be supplied separately. The Docker example includes these executables. FIMO and `fasta-get-markov` are MEME Suite tools and must be available in the active environment or configured by path. Check the installed tools before running TF analyses.

The application sources `viz_obj_helpers.R` from the application directory. Make sure this helper file is present alongside `app.R`, because it provides the functions used to locate and read the per gene Feature Importance and epigenetic signal files.

### Prepare the included chromosome 21 example

The repository contains the original IHEC CRE-RF chromosome 21 object in [`epivoi_object/data/`](epivoi_object/data/) and, in the development branch, the 346 corresponding gene-specific Feature Importance files and 346 epigenetic signal files. The original RDS still contains paths from the research server. After cloning, run this command **from the repository root**:

```bash
Rscript scripts/prepare_chr21_example.R
```

The script checks the supplied matrix files and creates `epivoi_object/data/epivoi_chr21_local.rds` with portable relative file references. **The original `epivoi_obj_cre_rf_IHEC_1MB_test_performance_chr21.rds` is not modified**, and the generated local RDS should not be committed.

See the [chromosome 21 example instructions](epivoi_object/data/README.md) for file layout and usage.

### Configure resource locations

Copy the provided configuration template, then edit the local copy:

```bash
cp .Renviron.example .Renviron
```

The example configuration below combines the **actual repository path of the local chromosome 21 object** with **actual paths used on the research server** for external resources. Replace the `/projects/...` values with paths to resources on your machine. The external files and executables are not automatically installed by cloning the repository.

```text
EPIVOI_OBJECT_PATH=epivoi_object/data/epivoi_chr21_local.rds
EPIVOI_REFSEQ_FILE=/projects/single_cell_stitchit/work/visualization/app/local_files/ncbiRefSeq.txt.gz
EPIVOI_IGV_ASSETS_DIR=/projects/single_cell_stitchit/work/visualization/app/igv
EPIVOI_TRAP_BIN=/projects/single_cell_stitchit/work/visualization/app/PASTAA_local/TRAP
EPIVOI_PASTAA_BIN=/projects/single_cell_stitchit/work/tf_enrichment/PASTAA/PASTAA
EPIVOI_PASTAA_ENERGY_MATRIX=/projects/single_cell_stitchit/work/visualization/app/tf_pastaa_output/matrices/Jaspar_Hocomoco_Kellis_human_energy.txt
EPIVOI_GENOME_FASTA=/projects/single_cell_stitchit/work/visualization/app/local_resources/hg38.fa
EPIVOI_MEME_MOTIF_FILE=/projects/single_cell_stitchit/work/visualization/tf_enrichment/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.meme
EPIVOI_FIMO_ENCODE_BG=/path/to/compatible/fimo_background.txt
EPIVOI_PASTAA_CORES=4
```

`EPIVOI_REGULATORY_BUILD_BED` is optional. Set it to the path of a compatible local regulatory annotation BED file if you want to use a local track rather than dynamic Ensembl Regulatory Build retrieval, for example `EPIVOI_REGULATORY_BUILD_BED=/path/to/regulatory_build.bed`.

The example values identify the original research configuration, not bundled input files. In particular, the correct FIMO background file must be verified for your run. A MEME motif database is **not** interchangeable with a standalone FIMO background frequency file even though MEME files can contain background frequency headers.

**IGV.js assets are external JavaScript and CSS inputs.** Obtain compatible `igv.min.js` and `igv.min.css` files and place them together in a local directory. Set `EPIVOI_IGV_ASSETS_DIR` to that directory. The code loads them as `igv/igv.min.js` and `igv/igv.min.css` through its Shiny resource mapping. If you store these resources elsewhere, update the configured local path, not the JavaScript URL in the browser.

Do not commit `.Renviron` containing machine specific paths or credentials. Values in this file are read when R starts, so restart the application after changes.

### Start EpIVoi with a new object

Start from the repository root. After preparing the chromosome 21 example, run:

```bash
Rscript app.R epivoi_object/data/epivoi_chr21_local.rds
```

To start with another compatible input object, pass its actual path instead:

```bash
Rscript app.R /absolute/path/to/your_epivoi_object.rds
```

Alternatively, set `EPIVOI_OBJECT_PATH` in `.Renviron` and run:

```bash
Rscript app.R
```

The application uses the supplied command line object path in preference to `EPIVOI_OBJECT_PATH`. Open the local address printed by Shiny, or access it through a suitable SSH tunnel when running on a remote server.

For detailed visual guidance, see the [module screenshots](docs/screenshots/). A separate screenshot-matched case study tutorial can be linked here when available. The main README intentionally does not repeat a long typical workflow.

## Input

EpIVoi expects outputs from a gene-wise prediction model and **does not train models**. The following section gathers **all input files and external resources in one place**, including the optional Regulatory Build annotation.

### Included chromosome 21 input dataset

The documented working example uses the original [IHEC CRE-RF chromosome 21 object](epivoi_object/data/epivoi_obj_cre_rf_IHEC_1MB_test_performance_chr21.rds). It covers **346 genes** on human **hg38**. The additional data are included in these repository directories:

| Input | Repository location | Quantity |
| --- | --- | ---: |
| Original EpIVoi object | [`epivoi_object/data/epivoi_obj_cre_rf_IHEC_1MB_test_performance_chr21.rds`](epivoi_object/data/epivoi_obj_cre_rf_IHEC_1MB_test_performance_chr21.rds) | 1 |
| Feature Importance matrices | [`epivoi_object/data/feature_importance/`](epivoi_object/data/feature_importance/) | 346 |
| Epigenetic signal matrices | [`epivoi_object/data/epigenetic_signal/`](epivoi_object/data/epigenetic_signal/) | 346 |
| Local object preparation | [`scripts/prepare_chr21_example.R`](scripts/prepare_chr21_example.R) | 1 script |

Run `Rscript scripts/prepare_chr21_example.R` once after cloning. The original object references absolute server paths, so the generated `epivoi_chr21_local.rds` is necessary for portable use of the bundled matrices. The preparation script preserves all 346 genes and does not replace the original object. See [example setup details](epivoi_object/data/README.md).

### Standardized EpIVoi object

The main input is a serialized R object (`.rds`) with S3 class `epivoi_object`. It is independent of the prediction algorithm and includes the following components:

```text
epivoi_object
    species
    genome
    gtf
    performance
        error_metric_name
        mse
        correlation
    feature_importance
        file_index
        groups
    epigenetic_signal
        file_index
        aggregation
```

`species` and `genome` identify the organism and assembly. `gtf` stores gene annotation. The `performance` component contains per-gene correlation and error metrics. Each `file_index` maps `gene_id` to the corresponding gene-specific data file. `feature_importance$groups` defines metadata groups such as cell types, conditions or their combinations. Only the requested per-gene files are loaded into memory during exploration.

The included `viz_obj_helpers.R` reads absolute paths as supplied and resolves relative paths against the directory of the loaded RDS when the application sets `object_dir`. Consequently, users can also supply their own compatible RDS with linked matrices. The chromosome 21 preparation script converts the server-specific paths into this portable layout.

### Per gene Feature Importance matrices

The current project Feature Importance files are tab separated (`.tsv`). Their first column identifies the biological group, while remaining columns contain genomic region coordinates and numerical feature importance values. In the current helper implementation, the matrix is transposed on reading so downstream modules work with regions as rows and groups as columns. Positive and negative values must be retained.

See the gene-specific files in [`epivoi_object/data/feature_importance/`](epivoi_object/data/feature_importance/) for real examples of the expected format.

### Per gene epigenetic signal matrices

The epigenetic signal data can be stored as compressed tab separated files (`.tsv.gz`). The first column is `region`, with additional columns for biological groups and numerical signal values. These signals describe the epigenetic features supplied to the original gene model, for example ATAC accessibility or H3K27ac levels.

See [`epivoi_object/data/epigenetic_signal/`](epivoi_object/data/epigenetic_signal/) for matching compressed input files.

### Annotations, reference sequences and motif files

The gene-wise matrices are bundled for chromosome 21, but the following resources remain **external** unless a compatible file has been explicitly added to the repository. Use the paths in `.Renviron.example` as a guide and replace them with paths valid on your machine.

| Resource | Environment variable | Expected content and example source path |
| --- | --- | --- |
| RefSeq gene annotations | `EPIVOI_REFSEQ_FILE` | Gzipped UCSC RefSeq annotation table, e.g. `/projects/single_cell_stitchit/work/visualization/app/local_files/ncbiRefSeq.txt.gz` |
| Reference genome | `EPIVOI_GENOME_FASTA` | Indexed genome FASTA consistent with the object's `genome` (`hg38` for chromosome 21), e.g. `/projects/single_cell_stitchit/work/visualization/app/local_resources/hg38.fa` |
| TRAP binding energy matrix | `EPIVOI_PASTAA_ENERGY_MATRIX` | TRAP-compatible energy matrix with parameter header and numerical energy values, e.g. `/projects/single_cell_stitchit/work/visualization/app/tf_pastaa_output/matrices/Jaspar_Hocomoco_Kellis_human_energy.txt` |
| MEME motif collection | `EPIVOI_MEME_MOTIF_FILE` | MEME-format TF motifs, e.g. `/projects/single_cell_stitchit/work/visualization/tf_enrichment/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.meme` |
| FIMO background | `EPIVOI_FIMO_ENCODE_BG` | MEME-format motif resource used by the current ENCODE-mode implementation |
| Ensembl Regulatory Build | `EPIVOI_REGULATORY_BUILD_BED` (optional) | Local BED annotation with chromosome, start and end, plus optional annotation fields |

The RefSeq and FASTA coordinate system must match the object assembly. **A MEME motif database is not interchangeable with the stand-alone FIMO background file.** When a local Regulatory Build BED is not specified, the app can retrieve compatible Ensembl annotations dynamically with `biomaRt`, which requires network access.

### Executables and browser assets

| Resource | Setting | Requirement |
| --- | --- | --- |
| IGV.js browser assets | `EPIVOI_IGV_ASSETS_DIR` | Directory containing `igv.min.js` and `igv.min.css`. Files are supplied externally and the configured directory must match their local storage location. |
| TRAP | `EPIVOI_TRAP_BIN` | Locally accessible TRAP executable. |
| PASTAA | `EPIVOI_PASTAA_BIN` | Locally accessible PASTAA executable. |
| FIMO | `EPIVOI_FIMO_BIN` (optional) | FIMO executable, otherwise searched on `PATH`. |
| fasta-get-markov | `EPIVOI_FASTA_GET_MARKOV_BIN` (optional) | MEME Suite executable, otherwise searched on `PATH`. |

The chromosome 21 dataset supplies the RDS and per-gene matrices but **does not bundle TRAP, PASTAA, IGV.js, complete hg38 or the external motif resources**. Those must be installed or supplied separately and configured before enabling their modules.

## Functions and modules

The functional grouping follows the EpIVoi manuscript: performance assessment, Feature Importance, locus exploration in IGV and TF enrichment. The names below reflect the application's current interface.

### Performance Overview

This module summarizes prediction quality across gene models using the available error metric and Pearson or Spearman test correlations. The distribution plot can display a density curve or histogram. Correlation based classes separate failed, low, medium and high performing models. Performance data are searchable and downloadable.

**Parameters:** **Metric for distribution plot** selects the error or test correlation measure. **Distribution type** chooses histogram or density. **Correlation metric for class cutoffs** selects the metric used to classify genes. **Low/medium cutoff** and **Medium/high cutoff** set the two positive correlation boundaries. **Show correlation classes** filters displayed classes. **Search gene** filters by gene symbol or identifier. Models with missing or nonpositive test correlation fall into the failed class.

[Overview screenshot](docs/screenshots/01_overview.pdf)

#### Train vs Test

The dedicated Train vs Test view compares training and test performance where both measurements are available. It can display paired errors or supported correlations and highlight the best performing genes.

**Parameters:** **Metric for train vs test scatter** selects the paired metric. **Use log10 scale** transforms the error display when appropriate. **Show test correlation classes** filters gene classes. **Correlation metric for top gene highlighting** and **Number of top genes to highlight** determine highlighted genes. Tables can be downloaded.

[Train vs Test screenshot](docs/screenshots/02_train_test.pdf)

### Feature Importance

The Feature Importance module displays regions contributing to prediction of a chosen gene within a chosen biological group. Ranking uses the magnitude of Feature Importance, while the original sign is kept in plots and tables. Feature Importance is model dependent and does not by itself establish direct regulation.

**Parameters:** **Biological group** selects the metadata-defined biological group encoded in the object. **Gene symbol** selects the gene. **Show top N regions by absolute Feature Importance** sets the number of regions to display. **Filter regions by minimum absolute Feature Importance** removes lower magnitude regions.

[Feature Importance screenshot](docs/screenshots/03_feature_importance.pdf)

### Genome Browser (IGV)

The embedded IGV genome browser displays gene annotations, Feature Importance and epigenetic signal at selected loci. Tracks can be compared across biological groups. Optional Ensembl Regulatory Build annotations and FIMO predicted TF binding sites support genomic interpretation. Selected tracks can be exported as vector PDF.

**Parameters:** **Load** updates the browser. **Gene** identifies the locus. **Cell types** selects biological groups for tracks. **Feature Importance type** and **Epigenetic Signal type** set user facing track labels. **Back to previous gene** navigates session history. **Show Ensembl regulatory build track** controls regulatory annotations. The FIMO controls are described immediately below, because FIMO is launched from this module.

#### FIMO motif analysis in IGV

FIMO scans candidate genomic sequences for occurrences of selected TF motifs and shows the resulting genomic positions as tracks. It uses the configured MEME motif collection, reference genome and nucleotide background resource.

**User controlled parameters in the current interface:**

| Parameter | Behavior |
| --- | --- |
| TF source for FIMO | Uses the latest saved PASTAA TFs or manual TF/motif names. |
| Manual TF names | Accepts up to 20 TF names or motif identifiers when manual mode is chosen. |
| Top saved PASTAA TFs to show in FIMO | Controls the number of saved TFs used in saved mode, from 1 to 20. |
| Replace previous FIMO tracks | Replaces previous predicted binding site tracks or retains them. |

**Fixed settings in the current implementation:** FIMO runs with `--motif-pseudo 0.1` and `--thresh 1e-4`. The threshold is a **p value threshold**, because `--qv-thresh` is not enabled. The interface currently does **not** offer a user parameter for choosing between p value and q value thresholds or for changing this cutoff. These values must be changed in the code if another threshold is needed. PASTAA's user selected **Max q-value / BH FDR** is a distinct control, not a FIMO q value cutoff.

[IGV screenshot](docs/screenshots/04_igv.pdf) · [FIMO screenshot](docs/screenshots/06_fimo.pdf)

### Transcription Factor Enrichment

EpIVoi integrates TRAP binding affinity calculations with PASTAA enrichment testing. Regions are ordered according to their Feature Importance and used to identify TFs associated with the prioritization. Direction reflects whether a region contributes positively or negatively to the model prediction and **must not** be interpreted as a TF's activating or repressing role.

**Run per gene** ranks regions for a chosen gene separately in each selected biological group. **Run per biological group** pools regions across eligible genes. To make values comparable, each gene's Feature Importance values are converted to within gene z scores before pooling. Duplicate genomic regions are removed. Filtering genes by test correlation before loading region files can reduce runtime. Significance is evaluated using PASTAA p values and Benjamini Hochberg adjusted q values.

**Parameters:** **PASTAA run mode** selects per gene or per biological group. **Gene** is used in gene mode. **Select biological groups** chooses one or more groups. **Exclude low-performing gene models**, **Performance metric** and **Minimum selected test correlation** control optional filtering in group mode. **Number of top regions** sets the size of each ranked input set. **Feature Importance direction** chooses positive, negative, both directions separately or absolute ranking. **Max q-value / BH FDR** sets the significance filter. **Save top N TFs** controls the number of TFs saved for downstream visualization in FIMO.

Results can be visualized as TF bar plots and comparisons of positive and negative groups, with table and figure exports. Absolute Feature Importance results are presented separately.

#### Parallel execution

In **biological group mode**, loading Feature Importance input files is parallelized **across genes** on Unix like systems. This input preparation is distinct from the later TRAP and PASTAA calculations.

In the current development version, prepared independent jobs for the selected biological groups and Feature Importance directions can also be processed in parallel. Preparation of the ranked regions is completed in the main Shiny process before executing parallel TRAP and PASTAA jobs. Results are collected afterward. The job execution therefore still blocks the current Shiny session while the batch completes.

Set the number of workers before launching the app:

```bash
export EPIVOI_PASTAA_CORES=4
```

Alternatively set `EPIVOI_PASTAA_CORES=4` in `.Renviron`. On supported Unix like systems, the application uses fork based parallel execution when multiple workers and independent tasks are available. Other systems use a sequential fallback. Runtime and memory usage depend on the number of genes, tasks and available server resources.

[TF enrichment screenshot](docs/screenshots/05_tf_enrichment.pdf)

### Reproducibility and Help

The **Reproducibility** tab generates downloadable R code reflecting the current analysis configuration. It documents Feature Importance loading, group based z score normalization, region ranking, BED and FASTA preparation, TRAP and PASTAA analysis and BH FDR adjustment. Paths and outputs in the downloaded script depend on the local working directory and configured resources.

The **Help** tab provides in app explanations of the plots, analysis modes and interpretation of results.

## Implementation and external tools

EpIVoi is implemented in R/Shiny. Locus visualization uses IGV.js with Shiny to JavaScript messaging. Regulatory annotations may be retrieved using biomaRt. TF analysis uses TRAP, PASTAA and FIMO from MEME Suite. The accompanying manuscript reports the implementation versions IGV.js 3.3.1 and FIMO/MEME Suite 5.5.7. Your locally installed software versions should be checked independently.

## Screenshots and further documentation

The repository screenshot directory is [`docs/screenshots/`](docs/screenshots/). Documentation links above correspond to the relevant module. A detailed case study tutorial, if published, should use the same screenshot sequence and the same application version. This avoids duplicating a long typical workflow in the main README.

## License

Consult the repository for the applicable license and distribution terms. Do not assume a license if none has been selected.

#### Gene Features in Overview

The Overview tab includes four violin plots comparing gene features
across correlation classes (failed, low, medium and high):

- Number of transcription start sites (TSSs)
- Total exon length
- Gene length
- Gene density (±500 kb around the TSS)

Gene features were calculated using
[FuFis](https://github.com/DennisHeck/FuFis),
specifically `GTF_Processing.gene_feature_table()`.

The plots follow the existing correlation thresholds and class filters,
support PDF downloads, and use the gene feature table associated with
the loaded EpIVoi object.

See [DOCKER_SETUP.md](DOCKER_SETUP.md) for configuration details.
