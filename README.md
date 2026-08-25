# EpIVoi

EpIVoi is an R/Shiny application for interactive exploration of
gene-level model performance, feature importance, epigenetic signal,
transcription-factor enrichment and predicted TF binding sites.

## Main functionality

- Model performance overview
- Train vs test performance
- Feature Importance visualization
- IGV genomic tracks
- Epigenetic signal tracks
- PASTAA TF enrichment
- FIMO TF binding-site prediction
- Reproducibility code generation

## Feature Importance

The EpIVoi object contains raw Feature Importance / SHAP values.

For gene-level TF enrichment, raw SHAP values are used directly.

For cell-type-level TF enrichment, raw SHAP values are standardized
within each gene across genomic bins/regions before regions are pooled
and ranked across genes:

```r
score_z = as.numeric(scale(score))
