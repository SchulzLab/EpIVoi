# Full EpIVoi runtime: R/Shiny and MEME Suite (FIMO, fasta-get-markov).
# Includes the bundled nine-gene example and approved TF binaries.
FROM mambaorg/micromamba:1.5.10

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# Prepare the included sparse example reference and its FASTA index.
RUN gzip -dc /app/example_input_files/example_genome_sparse.fa.gz \
      > /app/example_input_files/example_genome_sparse.fa && \
    /opt/conda/bin/samtools faidx \
      /app/example_input_files/example_genome_sparse.fa && \
    test -s /app/example_input_files/example_genome_sparse.fa.fai


RUN mkdir -p /app/shap_bw \
    /app/fimo \
    /app/tf_pastaa_app_files/tmp \
    /app/tf_pastaa_app_files/results \
    /app/tf_pastaa_app_files/fimo && \
    chown -R $MAMBA_USER:$MAMBA_USER \
    /app/shap_bw /app/fimo /app/tf_pastaa_app_files

ENV EPIVOI_FIMO_BIN=/opt/conda/bin/fimo \
    EPIVOI_FASTA_GET_MARKOV_BIN=/opt/conda/bin/fasta-get-markov \
    EPIVOI_TRAP_BIN=/app/docker_example/bin/TRAP \
    EPIVOI_PASTAA_BIN=/app/docker_example/bin/PASTAA \
    EPIVOI_OBJECT_PATH=/app/example_input_files/example_epivoi_object.rds \
    EPIVOI_REFSEQ_FILE=/app/example_input_files/example_ncbiRefSeq.txt.gz \
    EPIVOI_GENOME_FASTA=/app/example_input_files/example_genome_sparse.fa \
    EPIVOI_IGV_ASSETS_DIR=/app/docker_example/igv \
    EPIVOI_PASTAA_ENERGY_MATRIX=/app/example_input_files/example_energy_matrix.txt \
    EPIVOI_MEME_MOTIF_FILE=/app/example_input_files/example_motifs.meme \
    EPIVOI_FIMO_ENCODE_BG=/app/example_input_files/example_motifs.meme \
    EPIVOI_PASTAA_CORES=2



EXPOSE 3838
CMD ["Rscript", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
