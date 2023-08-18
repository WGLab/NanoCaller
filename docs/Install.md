# Installation
There are three ways to install and run NanoCaller, via Docker, Singularity or Conda.

NanoCaller has been developed and tested to work with Linux OS; we do not recommend using Windows or Mac OS. However, if you use Windows or Mac OS and have Docker installed on your machine, you can run NanoCaller inside a Docker container. NanoCaller does not require a GPU or any other special hardware to run.

Please check the [NanoCaller Docker Hub repository](https://hub.docker.com/repository/docker/genomicslab/nanocaller) for the most up to date version of NanoCaller docker image.


## Conda Installation

If you do not have Anaconda, you will need to install it first. Here, we show how to install Miniconda, a minimal installation of Anaconda, which is much smaller and has a faster installation:

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Go through all the prompts (installation in `$HOME` is recommended). The installation should take about 10 minutes, including the installation of Miniconda. 

### Bioconda
You can install NanoCaller in conda using Bioconda recipe:
`conda install -c bioconda nanocaller`

It is recommened that you install NanoCaller in a new conda environment to avoid any package conflict and use mamba for fast installation, in the following way:
```
conda create -n nanocaller_env -c conda-forge mamba
conda activate nanocaller_env
mamba install -c bioconda nanocaller
```

### Manual Installation
You can obtain the latest NanoCaller version from github that has not yet been pushed to bioconda via manual installation.

```
git clone https://github.com/WGLab/NanoCaller.git
conda env create -f NanoCaller/environment.yml
chmod +x NanoCaller/NanoCaller
conda activate nanocaller_env
```

Then you can run NanoCaller using `PATH_TO_NANOCALLER_REPO/NanoCaller --help`.


## Docker Installation
For instructions regarding Docker installation, please visit [Docker website](https://docs.docker.com/get-docker). There are three ways to obtain a Docker image for NanoCaller.

### 1) via Docker Hub (preferred)
You can pull NanoCaller docker images from Docker Hub by specifiying a version number.  
```
VERSION="3.4.1"
docker run genomicslab/nanocaller:${VERSION} NanoCaller --help
```

### 2) Locally build docker image
If you want to build an image for the latest commit of NanoCaller Github repository, use the following commands:

```
git clone https://github.com/WGLab/NanoCaller.git
docker build -t nanocaller NanoCaller
docker run nanocaller NanoCaller --help
```

## Singularity
For instructions regarding Singularity installation, please visit [Singularity website] (https://sylabs.io/guides/3.7/user-guide/quick_start.html).
```
VERSION="3.4.1"
singularity pull docker://genomicslab/nanocaller:${VERSION}
singularity exec -e --pwd /app nanocaller_${VERSION}.sif NanoCaller --help
```
