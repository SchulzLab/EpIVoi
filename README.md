# EpIVoi

EpIVoi is an interactive **R/Shiny application** for exploring model performance, Feature Importance, epigenetic signal, transcription factor enrichment, and genomic context from EpIVoi visualization objects.

The application combines model-level and gene-level summaries with interactive genomic visualization and downstream TF analysis.

[View full page EpIVoi Overview screenshot](https://github.com/SchulzLab/EpIVoi/blob/main/docs/screenshots/01_overview.pdf)

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

The environment contains R and the R packages required by the Shiny application together with MEME Suite utilities used by FIMO.

PASTAA and TRAP are external executables and must be available separately.

### 3. Configure external resources

Copy the provided environment template.

```bash
cp .Renviron.example .Renviron
```

Edit `.Renviron` and replace the placeholder paths with paths to the required external resources.

The local `.Renviron` file is ignored by Git and should not be committed.

Example configuration:

```text
EPIVOI_OBJECT_PATH=/path/to/epivoi_object.rds
EPIVOI_REFSEQ_FILE=/path/to/ncbiRefSeq.txt.gz
EPIVOI_IGV_ASSETS_DIR=/path/to/igv
EPIVOI_TRAP_BIN=/path/to/TRAP
EPIVOI_PASTAA_BIN=/path/to/PASTAA
EPIVOI_PASTAA_ENERGY_MATRIX=/path/to/energy_matrix.txt
EPIVOI_GENOME_FASTA=/path/to/reference_genome.fa
EPIVOI_MEME_MOTIF_FILE=/path/to/motifs.meme
EPIVOI_FIMO_ENCODE_BG=/path/to/fimo_motif_or_background_file
EPIVOI_REGULATORY_BUILD_BED=/path/to/regulatory_build.bed
EPIVOI_PASTAA_CORES=4
```

`EPIVOI_REGULATORY_BUILD_BED` is optional.

When it is not supplied, EpIVoi can retrieve compatible Ensembl Regulatory Build annotations dynamically for supported genomes.

`EPIVOI_OBJECT_PATH` is also optional when the object path is supplied directly at application startup.

### 4. Start EpIVoi with a new object

Provide the visualization object directly when starting the application.

```bash
Rscript app.R /path/to/epivoi_object.rds
```

Alternatively, set `EPIVOI_OBJECT_PATH` in `.Renviron` and start the application with:

```bash
Rscript app.R
```

A detailed walkthrough with matching screenshots is available in [`docs/toy_example.md`](docs/toy_example.md).

---

## Input

EpIVoi operates on outputs from gene-specific prediction models and does not perform model training.

### EpIVoi visualization object

The main application input is an RDS visualization object.

The object contains or references the data required by the enabled modules.

The object can contain:

* Species information
* Genome assembly information
* Gene annotation in GTF format
* Gene-level model performance
* Region-level Feature Importance
* Epigenetic signal
* File indexes and metadata used to locate data stored on disk

Feature Importance must be available for the biological groups that should be explored in the application.

A biological group can represent a cell type, disease condition, treatment, developmental stage, or a combination of metadata variables.

A minimal conceptual structure is:

```text
epivoi_object.rds
├── species
├── genome
├── gtf
├── performance
│   ├── mse
│   └── correlation
├── Feature Importance references
├── epigenetic signal references
└── file indexes
```

The exact object structure and helper functions used to construct an EpIVoi object are described in the detailed documentation and tutorial.

---

## External input files and resources

### RefSeq annotation

Environment variable:

```text
EPIVOI_REFSEQ_FILE=/path/to/ncbiRefSeq.txt.gz
```

Expected format: gzipped tab-separated RefSeq annotation table.

Example:

```text
585 NM_001001 chr1 + 11873 14409 11873 11873 3 11873,12612,13220, 12227,12721,14409, 0 DDX11L1
```

### IGV.js assets

Environment variable:

```text
EPIVOI_IGV_ASSETS_DIR=/path/to/igv
```

The directory must contain the JavaScript and CSS resources used by the embedded IGV browser.

Example:

```text
igv.min.js
igv.min.css
```

### TRAP executable

Environment variable:

```text
EPIVOI_TRAP_BIN=/path/to/TRAP
```

Expected input: executable TRAP binary.

### PASTAA executable

Environment variable:

```text
EPIVOI_PASTAA_BIN=/path/to/PASTAA
```

Expected input: executable PASTAA binary.

### TRAP binding energy matrix

Environment variable:

```text
EPIVOI_PASTAA_ENERGY_MATRIX=/path/to/energy_matrix.txt
```

Expected format: TRAP-compatible energy matrix containing motif definitions and the nucleotide background used for affinity calculation.

Example structure:

```text
>TF_NAME
A ...
C ...
G ...
T ...
```

### Reference genome

Environment variable:

```text
EPIVOI_GENOME_FASTA=/path/to/reference_genome.fa
```

Expected format: FASTA matching the genome assembly of the EpIVoi object.

Example:

```text
>chr1
NNNNACGTACGTACGT...
```

### FIMO motif collection

Environment variable:

```text
EPIVOI_MEME_MOTIF_FILE=/path/to/motifs.meme
```

Expected format: MEME motif file.

Example structure:

```text
MEME version 4

ALPHABET= ACGT

MOTIF TF_NAME

letter-probability matrix: alength= 4 w= ...
...
```

### FIMO motif and background resource

Environment variable:

```text
EPIVOI_FIMO_ENCODE_BG=/path/to/fimo_motif_or_background_file
```

The configured resource must be compatible with the FIMO workflow used by the application.

In the current submission workflow, the motif file contains the nucleotide background frequencies used by FIMO.

Example background header:

```text
Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25
```

### Regulatory Build

Optional environment variable:

```text
EPIVOI_REGULATORY_BUILD_BED=/path/to/regulatory_build.bed
```

Expected format: BED-compatible regulatory annotation with chromosome, start, and end coordinates.

Additional columns can contain names or annotation types.

Example:

```text
chr1    100000    100500    enhancer
chr1    105000    105300    promoter
```

When no local Regulatory Build file is supplied, EpIVoi attempts to retrieve compatible Ensembl regulatory annotations dynamically for supported genomes.

---

## Functions and modules

EpIVoi is organized into modules that correspond to the main analysis components described in the manuscript.

### 1. Overview

The **Overview** module summarizes model performance across all genes.

Users can inspect the distribution of a selected test performance metric, switch between density and histogram views, classify genes into performance classes, search for individual genes, and download the corresponding tables and plots.

#### Parameters

* **Metric for distribution plot** selects the displayed performance metric.
* **Distribution type** selects density or histogram visualization.
* **Correlation metric for class cutoffs** selects the test correlation used for classification.
* **Low/medium cutoff** defines the boundary between low and medium performance.
* **Medium/high cutoff** defines the boundary between medium and high performance.
* **Show correlation classes** controls which classes are displayed.
* **Search gene** filters the performance table.

Models with missing or nonpositive selected test correlation are assigned to the failed class.

[View full page Overview screenshot](https://github.com/SchulzLab/EpIVoi/blob/main/docs/screenshots/01_overview.pdf)

---

### 2. Train vs Test

The **Train vs Test** module compares model performance on the training and test datasets across genes.

#### Parameters

* **Metric for train vs test scatter** selects the performance metric.
* **Use log10 scale** applies logarithmic display to the error metric when available.
* **Show test correlation classes** controls the displayed model classes.
* **Correlation metric for top gene highlighting** selects the ranking metric.
* **Number of top genes to highlight** controls how many genes are emphasized.

The corresponding gene table can be inspected and downloaded.

[View full page Train vs Test screenshot](https://github.com/SchulzLab/EpIVoi/blob/main/docs/screenshots/02_train_test.pdf)

---

### 3. Feature Importance

The **Feature Importance** module visualizes genomic regions that contribute most strongly to the prediction for a selected gene and biological group.

#### Parameters

* **Biological group** selects the metadata-defined group used for Feature Importance.
* **Gene** selects the gene to inspect.
* **Show top N regions by absolute Feature Importance** controls the number of displayed regions.
* **Minimum absolute Feature Importance** removes regions below the selected magnitude threshold.

Regions are ranked by the absolute magnitude of Feature Importance while the original sign is retained in the plot and table.

[View full page Feature Importance screenshot](https://github.com/SchulzLab/EpIVoi/blob/main/docs/screenshots/03_feature_importance.pdf)

---

### 4. IGV

The **IGV** module places the selected gene and model-derived regions into their genomic context.

The view can contain:

* Gene annotation
* Feature Importance tracks
* Epigenetic signal tracks
* Ensembl Regulatory Build annotations
* FIMO-predicted TF binding sites

#### Parameters

* **Gene** selects the genomic locus.
* **Biological groups** selects one or more metadata-defined groups.
* **Feature Importance type** labels the Feature Importance track.
* **Epigenetic signal type** labels the corresponding signal track.
* **Show Ensembl regulatory build track** adds regulatory annotations.
* **TF source for FIMO** selects saved PASTAA TFs or manually entered TF names.
* **Top saved PASTAA TFs to show in FIMO** controls how many saved TFs are scanned.
* **Replace previous FIMO tracks** controls whether existing FIMO tracks are replaced.

Selected IGV tracks can be exported as a vector PDF.

[View full page IGV screenshot](https://github.com/SchulzLab/EpIVoi/blob/main/docs/screenshots/04_igv.pdf)

---

### 5. TF enrichment

The **TF enrichment** module integrates TRAP and PASTAA to test transcription factor enrichment in regions ranked by Feature Importance.

Two analysis modes are available.

#### Run per gene

A selected gene is analyzed separately for one or more biological groups.

#### Run per biological group

Regions are pooled across eligible genes for each selected biological group.

Raw Feature Importance values are converted to within-gene z-scores before pooling so that genes with different Feature Importance scales can be compared.

Gene models can optionally be filtered using a selected test correlation metric.

#### Parameters

* **PASTAA run mode** selects gene mode or biological group mode.
* **Gene** selects the gene in gene mode.
* **Exclude low performing gene models** activates performance filtering in biological group mode.
* **Performance metric** selects the test correlation used for filtering.
* **Minimum selected test correlation** sets the model performance threshold.
* **Select biological groups** selects the groups to analyze.
* **Save top N TFs** controls how many enriched TFs are retained for downstream use.
* **Max q value / BH FDR** filters PASTAA results after Benjamini-Hochberg correction.
* **Feature Importance direction** selects positive, negative, both separately, or absolute Feature Importance.
* **Number of top regions** defines how many ranked regions are used in each PASTAA analysis.

Positive and negative TF enrichment can be compared in a combined dot plot.

Absolute Feature Importance results are displayed separately.

#### Parallel Feature Importance loading

Biological group TF enrichment can require loading Feature Importance files for many genes.

On Unix-like systems, EpIVoi parallelizes this input loading step across genes.

The number of worker processes can be configured with:

```text
EPIVOI_PASTAA_CORES=4
```

If no value is supplied, EpIVoi selects a worker count automatically.

A sequential fallback is used when fork-based parallel processing is unavailable.

[View full page TF enrichment screenshot](https://github.com/SchulzLab/EpIVoi/blob/main/docs/screenshots/05_tf_enrichment.pdf)

---

### 6. FIMO motif analysis

FIMO is launched from the **IGV** module and identifies candidate TF binding sites in selected genomic regions.

TFs can be obtained from the latest saved PASTAA results or entered manually.

EpIVoi maps TF names to available motifs and adds successful FIMO predictions to the IGV view.

#### User-controlled parameters

* **TF source for FIMO** selects saved PASTAA TFs or manual TF input.
* **Manual TF names** accepts up to 20 TF or motif names.
* **Top saved PASTAA TFs to show in FIMO** sets the number of saved PASTAA TFs used when the saved source is selected.
* **Replace previous FIMO tracks** determines whether previous predictions are replaced.

#### Current FIMO settings

The current implementation calls FIMO with:

* **Motif pseudocount:** `0.1`
* **Significance threshold:** `1e-4`
* **Threshold type:** P-value

The threshold is interpreted as a P-value threshold because `--qv-thresh` is not enabled in the current implementation.

The motif collection and nucleotide background are provided through the configured external resources.

[View full page FIMO screenshot](https://github.com/SchulzLab/EpIVoi/blob/main/docs/screenshots/06_fimo.pdf)

---

### 7. Reproducibility

The **Reproducibility** module generates an R script that reflects the current EpIVoi analysis settings.

The generated script documents the major analysis steps outside the GUI, including:

* Feature Importance loading
* Biological group z-score calculation
* Region ranking
* BED and FASTA preparation
* TRAP
* PASTAA
* BH FDR correction

The script itself is downloaded through the browser using **Download reproducibility code**.

When the generated script is executed, analysis files are written relative to the R working directory unless the user supplies a path in the output prefix.

Generated files include:

* BED files
* Ranked region files
* FASTA files
* TRAP affinity files
* PASTAA output files

---

### 8. Help

The **Help** module contains an in-app explanation of the main EpIVoi concepts, plots, analysis modes, and interpretation of results.

---

## Software versions used for the implementation

The current EpIVoi implementation was tested with the following external software:

* **IGV.js 3.3.1** for interactive genomic visualization
* **MEME Suite / FIMO 5.5.7** for prediction of transcription factor motif occurrences
* **TRAP** local executable for calculation of transcription factor binding affinities
* **PASTAA** local executable for transcription factor enrichment analysis

IGV is implemented directly using IGV.js and is not based on the `igvShiny` R package.

---

## Documentation screenshots

Documentation screenshots are stored in:

```text
docs/screenshots/
```

Expected files:

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

Add the project license here if and when a license is selected.

