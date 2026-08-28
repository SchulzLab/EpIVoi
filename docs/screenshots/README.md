# EpIVoi screenshot folder

Place the six final screenshots in this directory using exactly these filenames:

```text
01_overview.pdf
02_train_test.pdf
03_feature_importance.pdf
04_igv.pdf
05_tf_enrichment.pdf
06_fimo.pdf
```

The main `README.md` already references these paths, so GitHub will display the screenshots automatically once the PNG files are added here.

## Screenshot checklist

### `01_overview.pdf`

Show:

- the Overview tab;
- the navigation bar;
- representative performance plots;
- preferably a gene search or part of the performance table.

### `02_train_test.pdf`

Show:

- the Train vs Test tab;
- the train/test scatter plot;
- representative highlighted genes.

### `03_feature_importance.pdf`

Show:

- a selected gene;
- a selected biological group;
- the Feature Importance bar plot;
- enough of the table or controls to identify the workflow.

### `04_igv.pdf`

Show:

- a loaded gene locus;
- gene annotation;
- SHAP / Feature Importance track;
- epigenetic-signal track;
- Ensembl regulatory annotations if enabled.

### `05_tf_enrichment.pdf`

Prefer a screenshot showing:

- TF-enrichment controls;
- combined positive/negative TF-enrichment dot plot;
- readable biological-group and TF labels.

If the current bar plot and absolute dot plot fit cleanly, they may also be included.

### `06_fimo.pdf`

Show:

- a successfully loaded gene locus in IGV;
- FIMO motif-site tracks;
- enough genomic context to identify the result.

## Style

- Use the same browser width where possible.
- Crop unnecessary browser chrome.
- Do not show usernames, private filesystem paths, internal hostnames, or terminal windows.
- Use representative non-empty results.
- Make sure labels remain readable when rendered on GitHub.
