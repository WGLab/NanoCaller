# Installation

NanoCaller has been developed and tested to work with Linux OS; we do not recommend using Windows or Mac OS. NanoCaller does not require a GPU or any other special hardware to run.
There are two ways to install and run NanoCaller, via docker or conda.

## Docker Installation
There are three ways to obtain a Docker image for NanoCaller.

### 1) via Docker Hub (preferred)
You can pull NanoCaller docker images from Docker Hub by specifiying a version number.  
```
VERSION="0.3.0"
docker run umahsn/nanocaller:${VERSION} python NanoCaller.py --help
```

### 2) Locally build docker image
If you want to build an image for the latest commit of NanoCaller Github repository, use the following commands:

```
git clone https://github.com/WGLab/NanoCaller.git
docker build -t nanocaller NanoCaller
docker run  nanocaller python NanoCaller.py --help
```

### 3) Saved docker image via Github release
If you want to use NanoCaller docker image saved in a tar file, download the image file by specifying a version number and use `docker load`.

```
VERSION="0.3.0"
wget https://github.com/WGLab/NanoCaller/releases/download/v${VERSION}/nanocaller_docker.tar.gz
docker load --input nanocaller_docker.tar.gz
docker run  nanocaller python NanoCaller.py --help
```


## Conda Installation
First, install Miniconda, a minimal installation of Anaconda, which is much smaller and has a faster installation.
Note that this version is meant for Linux below, macOS and Windows have a different script:

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Go through all the prompts (installation in `$HOME` is recommended). After Anaconda is installed successfully, simply run:

```
git clone https://github.com/WGLab/NanoCaller.git
cd NanoCaller
conda env create -f environment.yml
conda activate NanoCaller
```
The installation should take about 10 minutes, including the installation of Miniconda.
