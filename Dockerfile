FROM continuumio/miniconda3
USER root

WORKDIR /app

SHELL [ "/bin/bash", "--login", "-c" ]
COPY environment.yml .

RUN conda env create -f environment.yml && conda clean --yes --all

RUN conda init bash

RUN echo "conda activate nanocaller_env" > ~/.bashrc

ENV PATH=/opt/conda/envs/nanocaller_env/bin:$PATH

COPY ./nanocaller_src /opt/conda/envs/nanocaller_env/bin/nanocaller_src
COPY ./NanoCaller /opt/conda/envs/nanocaller_env/bin
COPY ./NanoCaller_WGS /opt/conda/envs/nanocaller_env/bin

RUN chmod +x /opt/conda/envs/nanocaller_env/bin/NanoCaller
