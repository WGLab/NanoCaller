# Installation

NanoCaller has been developed and tested to work with Linux OS; we do not recommend using Windows or Mac OS. NanoCaller does not require a GPU or any other special hardware to run.
There are two ways to install and run NanoCaller, via docker or conda.
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
