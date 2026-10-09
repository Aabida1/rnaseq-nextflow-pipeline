FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER env/environment.yml /tmp/environment.yml
RUN micromamba create --yes --name pipeline --file /tmp/environment.yml && micromamba clean --all --yes
ENV PATH="/opt/conda/envs/pipeline/bin:$PATH"
WORKDIR /work
