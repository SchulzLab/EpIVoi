# EpIVoi — Example Input Files

This directory contains small examples illustrating the input file formats used by EpIVoi.

The examples are intended for documentation and format inspection. They are not complete datasets for biological analysis or a standalone runnable installation.

All example scores, signals, motif values, and genomic intervals are illustrative or synthetic unless otherwise indicated. The RefSeq example contains a small excerpt of reference annotation.

| File | Description |
|---|---|
| `example_epivoi_object.rds` | Example R object with performance metrics, GTF annotations, and SHAP/ATAC file indexes |
| `example_feature_importance.tsv` | Gene-level Feature Importance matrix, with biological groups as rows and genomic regions as columns |
| `example_epigenetic_signal.tsv.gz` | Epigenetic signal matrix, with genomic regions as rows and biological groups as columns |
| `example_ncbiRefSeq.txt.gz` | Example RefSeq gene annotation records |
| `example_genome.fa` | Short synthetic reference sequence in FASTA format |
| `example_energy_matrix.txt` | Illustrative PASTAA/TRAP energy matrix structure |
| `example_motifs.meme` | Synthetic transcription factor motif in MEME format |
| `example_fimo_background.txt` | Example zero-order nucleotide background frequencies for FIMO |
| `example_regulatory_build.bed` | Optional regulatory annotation example in BED format |

The example RDS references the accompanying SHAP and ATAC example files using relative paths. Applications loading the object must resolve these paths relative to the RDS location.

The energy matrix and reference files are format illustrations only and have not been validated as substitutes for complete analysis resources.

For real analyses, users must provide appropriate reference files, motif databases, executables, and EpIVoi visualization objects as described in the main repository README.
