# EpIVoi toy example

This walkthrough demonstrates a typical EpIVoi analysis using a single gene and one biological group.

The exact genes and biological groups available in the interface depend on the loaded visualization object. The example below uses:

```text
Gene: ABHD5
Biological group: Astrocytes:Alzheimer's
```

when these are available in the loaded object.

---

## 1. Start EpIVoi

Activate the Conda environment:

```bash
conda activate epivoi
```

Configure the required external resources and start the application:

```bash
Rscript app.R
```

Open the Shiny URL printed in the terminal.

---

## 2. Inspect the Overview tab

Open **Overview**.

The Overview tab provides a gene-level summary of model performance.

Try the following:

1. Select the desired performance metric.
2. Switch between the density and histogram views.
3. Inspect the violin plot and correlation-class bar plot.
4. Use the default correlation-class cutoffs or modify them.
5. Search for:

```text
ABHD5
```

The gene search accepts a gene symbol and, where available, an Ensembl ID.

The performance table can be downloaded for downstream analysis.

Screenshot reference:

```text
docs/screenshots/01_overview.png
```

---

## 3. Compare Train vs Test performance

Open **Train vs Test**.

This tab compares training and test performance across genes.

1. Select a train/test metric.
2. Inspect the scatter plot.
3. Optionally change the displayed correlation classes.
4. Choose a correlation metric for top-gene highlighting.
5. Increase or decrease the number of highlighted genes.

This view is useful for identifying genes whose test performance is strong enough for detailed downstream inspection.

Screenshot reference:

```text
docs/screenshots/02_train_test.png
```

---

## 4. Inspect Feature Importance

Open **Feature Importance**.

Select:

```text
Gene: ABHD5
Biological group: Astrocytes:Alzheimer's
```

Then:

1. Choose the number of top regions to display.
2. Optionally set a minimum absolute Feature Importance threshold.
3. Inspect the bar plot.
4. Inspect the Feature Importance table.

Regions are ranked by absolute Feature Importance, while the original signed score is retained.

Screenshot reference:

```text
docs/screenshots/03_feature_importance.png
```

---

## 5. Load the gene in IGV

Open **IGV**.

1. Select or enter `ABHD5`.
2. Select `Astrocytes:Alzheimer's` or another available biological group.
3. Click **Load**.
4. Optionally enable **Show Ensembl regulatory build track**.

The IGV view may display:

- gene annotation;
- SHAP / Feature Importance tracks;
- epigenetic signal such as ATAC;
- Ensembl regulatory regions;
- FIMO tracks after motif scanning.

The regulatory-build track does not require a local BED file in the standard configuration. If no local override is supplied, EpIVoi can obtain compatible Ensembl regulatory annotations dynamically.

Screenshot reference:

```text
docs/screenshots/04_igv.png
```

---

## 6. Run PASTAA in gene mode

Open **TF enrichment**.

Choose:

```text
PASTAA run mode: Run per gene
Gene: ABHD5
Biological group: Astrocytes:Alzheimer's
```

For a compact example, use:

```text
Number of top regions: 50
```

Run positive Feature Importance first.

A successful terminal message resembles:

```text
Running PASTAA: ABHD5 | Astrocytes:Alzheimer's | pos | mode=gene
```

Then run negative Feature Importance.

A corresponding message resembles:

```text
Running PASTAA: ABHD5 | Astrocytes:Alzheimer's | neg | mode=gene
```

The current-selection bar plot should reflect the selected result direction.

---

## 7. Inspect the combined positive/negative TF-enrichment plot

After positive and negative PASTAA results exist, inspect:

```text
Positive vs. Negative Feature Importance (combined, by group)
```

Positive and negative results are shown together in a combined dot plot.

The plot allows the TF enrichment associated with both Feature Importance directions to be compared directly across the shown biological groups.

Screenshot reference:

```text
docs/screenshots/05_tf_enrichment.png
```

---

## 8. Run absolute Feature Importance analysis

Select:

```text
Feature Importance direction: Absolute
```

and start PASTAA again.

A successful run should contain:

```text
| abs | mode=gene
```

in the terminal output.

After an absolute result exists, the section:

```text
Absolute Feature Importance (by group)
```

shows the corresponding dot plot.

If no absolute run has been performed yet, the application reports:

```text
No absolute-direction TF enrichment results are available.
```

This is expected until an `abs` result exists.

---

## 9. Run PASTAA in biological-group mode

Switch to:

```text
PASTAA run mode: Run per biological group
```

No gene is selected in this mode.

Choose one or more biological groups.

Optionally enable:

```text
Exclude low-performing gene models
```

and use, for example:

```text
Performance metric: test_Pearson
Minimum selected test correlation: 0.2
```

In biological-group mode, raw Feature Importance values are standardized within gene before regions are pooled across genes.

This avoids directly pooling genes with very different Feature Importance scales.

On Unix-like systems, gene-level Feature Importance files can be read in parallel.

To explicitly limit the number of workers before starting the application:

```bash
export EPIVOI_PASTAA_CORES=4
```

---

## 10. Run FIMO and display motif sites in IGV

Before running FIMO, make sure the desired gene has already been loaded in IGV.

In the IGV sidebar, choose a TF source:

```text
Use latest saved PASTAA TFs
```

or:

```text
Type TF / motif IDs manually
```

If saved PASTAA TFs are used, select the number of top TFs.

Then click:

```text
Run FIMO / show binding sites in IGV
```

EpIVoi maps TF names to available motifs and runs FIMO using the configured motif database and background.

Successful motif sites are added to the IGV visualization.

Screenshot reference:

```text
docs/screenshots/06_fimo.png
```

---

## 11. Reproduce the analysis outside the GUI

Open **Reproducibility**.

The tab displays R code corresponding to the current application settings and major processing steps.

The generated code includes the logic for:

- loading the EpIVoi object;
- loading Feature Importance data;
- biological-group within-gene z-score calculation;
- region ranking;
- BED/FASTA preparation;
- TRAP;
- PASTAA;
- BH-FDR correction.

Use **Download reproducibility code** to save the generated script.

---

## 12. Download results

Depending on the active view, EpIVoi provides downloads for:

- model-performance tables;
- Feature Importance regions;
- plots;
- PASTAA result tables;
- the current TF-enrichment bar plot;
- the combined positive/negative TF-enrichment plot;
- the absolute TF-enrichment dot plot;
- IGV vector PDF export;
- reproducibility code.

---

## Expected result

A correctly configured installation should allow the user to:

- open all available EpIVoi tabs;
- search for a gene;
- inspect model performance;
- compare train and test performance;
- inspect Feature Importance regions;
- load the selected locus in IGV;
- display available genomic tracks;
- run positive and negative PASTAA analyses;
- compare positive and negative TF enrichment in a combined dot plot;
- obtain an absolute TF-enrichment plot after an absolute run;
- run biological-group TF enrichment;
- run FIMO when the required external resources are available;
- download results and reproducibility information.
