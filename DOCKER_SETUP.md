# EpIVoi Docker Setup

## Quick Start

EpIVoi includes a ready-to-run Docker example with nine gene models,
two biological groups, reference annotations, and TF analysis resources.

The example requires no manual preparation of input files.

### Requirements

- Docker Engine with Docker Compose, or Docker Desktop
- Support for Linux/amd64 containers

### Run the Example

Clone the repository and switch to the development branch:

```bash
git clone https://github.com/SchulzLab/EpIVoi.git
cd EpIVoi
git checkout aleksandra-epivoi-updates
Build and start the application:
docker compose up --build


Open http://localhost:3838 in your browser.

To stop the application:
docker compose down


The Docker build automatically extracts the bundled sparse reference
FASTA and creates its index.
Included Example
All example input data are located in example_input_files/.

The dataset includes nine gene models and two biological groups:
- Astrocytes:Alzheimer's
- Astrocytes:Unaffected
The example provides gene-specific feature importance and epigenetic
signal matrices, RefSeq annotations, a sparse hg38 reference, a TRAP
energy matrix, and a small MEME motif collection.

TRAP, PASTAA, and IGV.js assets are included in docker_example/.
FIMO and other required runtime dependencies are installed during
the Docker build.

The sparse reference preserves the original hg38 sequences at selected
genomic positions, while other positions are filled with N. It is not
a replacement for the complete hg38 reference genome.

The example motif collection is intentionally small. PASTAA results
may not pass the default BH FDR threshold. To explore the PASTAA-to-FIMO
workflow with this demonstration dataset, a less restrictive cutoff
may be selected before running PASTAA.

The example is intended to demonstrate application functionality,
not to establish biological findings.
Using Your Own Data
The default Docker configuration runs the bundled example.

For another dataset, provide a compatible EpIVoi RDS object together
with its referenced feature importance and epigenetic signal files,
genome reference, annotations, and TF resources.

Mount these files into the container using a Docker Compose override
and configure the corresponding EPIVOI_* environment variables.

Paths recorded in the RDS file indices must be accessible from inside
the container.

For Gene Features in Overview, provide a matching precomputed CSV next
to the RDS using the .gene_features.csv sidecar naming convention,
or configure EPIVOI_GENE_FEATURES_CSV.

The bundled nine-gene example does not currently include a matching
Gene Features CSV. The remaining Overview functionality is available.

See README.md for the object format and module descriptions.
Platform and Testing
The Docker application was built and launched successfully on an
Apple Silicon Mac using Linux/amd64 emulation.

IGV and the TRAP/PASTAA workflow were tested with the included
example. Loading saved PASTAA TFs for downstream FIMO analysis
was also tested.

Deployment on the laboratory servers has not been tested.
