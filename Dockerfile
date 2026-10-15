# Full EpIVoi runtime: R/Shiny and MEME Suite (FIMO, fasta-get-markov).
# TRAP/PASTAA are supplied as user-provided Linux amd64 binaries via bind mounts.
FROM mambaorg/micromamba:1.5.10

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

RUN mkdir -p /app/shap_bw \
    /app/fimo \
    /app/tf_pastaa_app_files/tmp \
    /app/tf_pastaa_app_files/results \
    /app/tf_pastaa_app_files/fimo && \
    chown -R $MAMBA_USER:$MAMBA_USER \
    /app/shap_bw /app/fimo /app/tf_pastaa_app_files

ENV EPIVOI_FIMO_BIN=/opt/conda/bin/fimo \
    EPIVOI_FASTA_GET_MARKOV_BIN=/opt/conda/bin/fasta-get-markov \
    EPIVOI_TRAP_BIN=/opt/epivoi/bin/TRAP \
    EPIVOI_PASTAA_BIN=/opt/epivoi/bin/PASTAA \
    EPIVOI_OBJECT_PATH=/data/object.rds \
    EPIVOI_REFSEQ_FILE=/resources/ncbiRefSeq.txt.gz \
    EPIVOI_GENOME_FASTA=/resources/hg38.fa \
    EPIVOI_IGV_ASSETS_DIR=/resources/igv \
    EPIVOI_PASTAA_ENERGY_MATRIX=/resources/energy_matrix.txt \
    EPIVOI_MEME_MOTIF_FILE=/resources/motifs.meme \
    EPIVOI_PASTAA_CORES=2

ENV PATH="/opt/conda/bin:${PATH}"

EXPOSE 3838
CMD ["Rscript", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
