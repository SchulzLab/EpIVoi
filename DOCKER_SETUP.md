# Docker deployment (full EpIVoi application)

This configuration packages the full Shiny application and the MEME Suite utilities (`fimo`, `fasta-get-markov`). It does **not** bundle the lab-specific TRAP/PASTAA executables, private reference files or user datasets.

## Prerequisites

- Docker Engine + Compose plugin, or Docker Desktop.
- Linux/amd64 support (Apple Silicon uses emulation).
- Permission to use and distribute the chosen TF binaries and data.

## Prepare inputs

Create `docker_data`, `docker_resources`, `docker_bin` beside `compose.yaml`.

- Put a **compatible EpIVoi RDS** at `docker_data/object.rds`, or point `EPIVOI_OBJECT_PATH` to another path **inside the container**, for example `/data/myobject.rds`.
- Put each referenced SHAP and epigenetic signal file into the mounted data tree. Ensure the RDS `file_index$file_path` entries resolve under the runtime `object_dir` (for relative paths) or as valid absolute container paths. Paths from the host (such as `/projects/...`) are not usable unless explicitly bind-mounted to the same location.
- Put `hg38.fa`, `ncbiRefSeq.txt.gz`, `energy_matrix.txt` and `motifs.meme` in `docker_resources` (or change the environment mappings). The optional `background.txt` is required only for the corresponding FIMO ENCODE configuration, which has not yet been validated.
- Put your licensed `TRAP` and `PASTAA` Linux/amd64 executables in `docker_bin` and mark them executable (`chmod +x docker_bin/{TRAP,PASTAA}`).
- Put IGV assets in `docker_resources/igv` with filenames expected by `app.R`.
- Optionally add regulatory build data and set `EPIVOI_REGULATORY_BUILD_BED` in Compose.

Files required for specific modules can be omitted if those modules are intentionally unused **only if** `app.R` permits startup without them.

## Start

```sh
docker compose build
docker compose up
```

Then visit http://localhost:3838 on the machine running Docker. Stop using `docker compose down`.

The default Compose configuration binds the application to `127.0.0.1:3838`. When Docker runs on a remote server, use an SSH tunnel or another approved access method. To use another local port, change the left-hand port in `compose.yaml`, for example `127.0.0.1:3839:3838`.

### FASTA index

Before running PASTAA or FIMO, ensure that `docker_resources/hg38.fa.fai` exists and matches the reference genome. The `/resources` volume is mounted read-only, so bedtools cannot create this index inside the running application container.

If the index is missing, it can be generated using the built Docker image:

```sh
docker compose build
docker compose run --rm --no-deps --user root \
  --entrypoint /opt/conda/bin/bedtools \
  -v "$PWD/docker_resources:/index_resources" \
  epivoi getfasta \
  -fi /index_resources/hg38.fa \
  -bed /dev/null \
  -fo /tmp/index_test.fa
Verify that docker_resources/hg38.fa.fai was created. Then start the application with docker compose up -d.
Functional validation
Open EpIVoi, load a compatible gene and biological group in IGV, run PASTAA, then use the saved PASTAA TFs to run FIMO. The full pipeline was successfully tested on an Apple Silicon Mac using Linux/amd64 emulation. The ENCODE-specific FIMO configuration has not been validated.

## Verify inside the container

```sh
docker compose run --rm epivoi Rscript -e 'library(shiny); library(data.table); cat("R runtime OK\n")'
docker compose run --rm epivoi fimo --version
docker compose run --rm epivoi fasta-get-markov --version
```

## Limitations

- The Docker image and the IGV, TRAP, PASTAA and FIMO workflow were tested successfully on an Apple Silicon Mac using Linux/amd64 emulation. Deployment on the lab servers has not yet been tested.
- A container is only fully portable if input paths and referenced files are accessible inside it. The example chromosome 21 preprocessing script can help adapt its own data, but does not automatically adapt other RDS objects.
- Large genomes and model data stay on host mounts and are not included in a publicly distributed image.

## Running Docker with the example dataset

The repository includes a small example dataset with nine genes and
two biological groups. This allows testing without the full model object.

### 1. Prepare the sparse reference

From the repository root:

    gzip -dkf example_input_files/example_genome_sparse.fa.gz
    samtools faidx example_input_files/example_genome_sparse.fa

The uncompressed reference occupies approximately 540 MB.

### 2. Prepare the Docker inputs

Run from the repository root:

    mkdir -p docker_data docker_resources docker_bin

    cp example_input_files/example_epivoi_object.rds docker_data/object.rds
    cp -r example_input_files/feature_importance docker_data/
    cp -r example_input_files/epigenetic_signal docker_data/

    cp example_input_files/example_ncbiRefSeq.txt.gz docker_resources/ncbiRefSeq.txt.gz
    cp example_input_files/example_genome_sparse.fa docker_resources/hg38.fa
    cp example_input_files/example_genome_sparse.fa.fai docker_resources/hg38.fa.fai
    cp example_input_files/example_energy_matrix.txt docker_resources/energy_matrix.txt
    cp example_input_files/example_motifs.meme docker_resources/motifs.meme

Provide the lab-specific TRAP and PASTAA Linux/amd64 executables:

    cp /path/to/TRAP docker_bin/TRAP
    cp /path/to/PASTAA docker_bin/PASTAA
    chmod +x docker_bin/TRAP docker_bin/PASTAA

Replace /path/to/ with the actual binary locations.

Also provide the required IGV.js assets in docker_resources/igv/.

### 3. Build and start

    docker compose build
    docker compose up -d

Open http://localhost:3838 on the Docker host.

For a remote server, use SSH port forwarding or another
approved access method.

### 4. Validate

Load a gene from the example dataset, inspect its feature importance
and epigenetic tracks, then test PASTAA and FIMO.

The example motif file contains only RUNX1. Therefore, enrichment
may finish without statistically significant results.

The full Docker workflow was tested on Apple Silicon using
Linux/amd64 emulation. Server deployment remains to be validated.
