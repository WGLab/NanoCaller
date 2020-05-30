# NanoCaller
NanoCaller is a computational method that integrates long reads in deep convolutional neural network for the detection of SNPs/indels from long-read sequencing data. NanoCaller uses long-range haplotype structure to generate predictions for each SNP candidate variant site by considering pileup information of other candidate sites sharing reads. Subsequently, it performs read phasing, and carries out local realignment of each set of phased reads and the set of all reads for each indel candidate variant site to generate indel calling, and then creates consensus sequences for indel sequence prediction.

## Citing NanoCaller
Please cite: Ahsan, Umair and Liu, Qian and Wang, Kai. NanoCaller for accurate detection of SNPs and small indels from long-read sequencing by deep neural networks. bioRxiv 2019.12.29.890418; doi: https://doi.org/10.1101/2019.12.29.890418

## Installation
### Method 1
First, install Miniconda, a minimal installation of Anaconda, which is much smaller and has a faster installation.
Note that this version is meant for Linux below, macOS and Windows have a different script:

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Go through all the prompts (installation in `$HOME` is recommended).
After Anaconda is installed successfully, simply run:

```
git clone https://github.com/WGLab/NanoCaller.git
cd NanoCaller
conda env create -f environment.yml
conda activate NanoCaller
```
### Method 2
If you do not wish to use Anaconda, simply install the packages in the file `environment.yml` using `pip`.

```
git clone https://github.com/WGLab/NanoCaller.git
cd NanoCaller
```

##Usage
