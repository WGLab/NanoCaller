FROM continuumio/miniconda3
USER root

WORKDIR /app

SHELL [ "/bin/bash", "--login", "-c" ]
COPY environment.yml .

RUN conda env create -f environment.yml && conda clean --yes --all

RUN conda init bash

RUN echo "conda activate NanoCaller" > ~/.bashrc

ENV PATH=/opt/conda/envs/NanoCaller/bin:$PATH

COPY ./scripts ./
