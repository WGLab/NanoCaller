# NanoCaller
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/nanocaller/README.html)

NanoCaller is a computational method that integrates long reads in deep convolutional neural network for the detection of SNPs/indels from long-read sequencing data. NanoCaller uses long-range haplotype structure to generate predictions for each SNP candidate variant site by considering pileup information of other candidate sites sharing reads. Subsequently, it performs read phasing, and carries out local realignment of each set of phased reads and the set of all reads for each indel candidate variant site to generate indel calling, and then creates consensus sequences for indel sequence prediction.

NanoCaller is distributed under the [MIT License by Wang Genomics Lab](https://wglab.mit-license.org/).

## Latest Updates
_**v3.3.0** (July 14 2023)_: Detailed description of SNP calls, including unfiltered SNP calls for variants determined to be false by NanoCaller, and inclusion of per-base probability output. Quality score has been adjusted to be on Phred scale.
_**v3.2.0** (May 14 2023)_: Support added for haploid variant calling which has significant improvement in recall for indel calling. New feature generation methods and models are are used for haploid SNP and indel calling. Now chrY and chrM are assumed to be haploid, with additional parameter --haploid_X to specify if chrX is haploid. Another parameter --haploid_genome can be used for haploid variant calling on all chromosomes.

_**v3.0.1** (March 14 2023)_ : Several critical bugs regarding coverage normalization and integer overflow fixed. These bug affected very low and high coverage sample. The normalization bug was only introduced in v3.0.0 so any samples processed before that should not have been affected. Whereas integer overflow bug was much older and it only was affecting sample with more than 256 coverage.

_**v3.0.0** (June 7 2022)_ : A major update in API with single entry point for running NanoCaller. Major changes in parallelization routine with GNU parallel no longer used for whole genome variant calling.

_**v2.0.0** (Feb 2 2022)_ : A major update in API and installation instructions, with release of bioconda recipe for NanoCaller. Added support for indel calling in case of poor or non-existent phasing.

_**v1.0.0** (Aug 8 2021)_ : First post-production release with citeable DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5176764.svg)](https://doi.org/10.5281/zenodo.5176764)

_**v0.4.1** (Aug 3 2021)_ : Fixed a bug causing slower runtime in whole genome variant calling mode. 

_**v0.4.0** (June 2 2021)_ : Added NanoCaller models trained on ONT reads basecalled with Guppy v4.2.2 and Bonito v0.30, as well as R10.3 reads. Added new NanoCaller models trained with long CCS reads (15-20kb library selection). Improved indel calling with rolling window for candidate selection which helps with indels in low complexity regions.


## Installation
NanoCaller can be installed using Docker or Conda. The easiest way to install is from the bioconda channel:

`conda install -c bioconda nanocaller`

or using Docker:

```
VERSION="3.3.0"
docker pull genomicslab/nanocaller:${VERSION}
```
Please refer to [Installation](docs/Install.md) for instructions regarding installing NanoCaller through other methods.

## Usage
General usage of NanoCaller is described in [Usage](docs/Usage.md). Some quick usage examples:

- `NanoCaller --bam YOUR_BAM --ref YOUR_REF --cpu 10` will run NanoCaller on whole genome using 10 parallel processes.
- `NanoCaller --bam YOUR_BAM --ref YOUR_REF --cpu 10 --mode snps` will only call SNPs.
- `NanoCaller --bam YOUR_BAM --ref YOUR_REF --cpu 10 --mode snps --phase` will only call SNPs and phase them, and will additionally phase the BAM file (under intermediate_phase_files subfolder split by chromosomes).
- `NanoCaller --bam YOUR_BAM --ref YOUR_REF --cpu 10 --haploid_genome` will run NanoCaller on whole genome under the assumption that the genome is haploid.
- `NanoCaller --bam YOUR_BAM --ref YOUR_REF --cpu 10 --regions chr22:20000000-21000000 chr21` will NanoCaller on chr21 and chr22:20000000-21000000 only.

For a comprehensive case study of variant calling on Nanopore reads, see [ONT Case Study](docs/ONT%20Case%20Study.md), where we describe end-to-end variant calling pipeline for using NanoCaller, where we start with aligning FASTQ files of HG002, calls variants using NanoCaller, and evaluate performances on various genomic regions.


## Trained models
Trained models for [ONT](https://github.com/WGLab/NanoCaller/tree/master/nanocaller_src/release_data/ONT_models) data, [CLR](https://github.com/WGLab/NanoCaller/tree/master/nanocaller_src/release_data/clr_models) data and [HIFI](https://github.com/WGLab/NanoCaller/tree/master/nanocaller_src/release_data/hifi_models) data can be found [here](https://github.com/WGLab/NanoCaller/tree/master/nanocaller_src/release_data). These models are trained on chr1-22 of the genomes stated below, unless mentioned othewise.

You can specify SNP and indel models using `--snp_model` and `--indel_model` parameters with a model name from tables below. For instance, if you want to use 'ONT-HG002\_bonito' SNP model and 'ONT-HG002' indel model, use the following command:

`NanoCaller --snp_model ONT-HG002_bonito --indel_model ONT-HG002`


#### SNP Models

| Model Name                 | Sequencing Technology | Genome          | Coverage | Benchmark | Basecaller   |
| -------------------------- | --------------------- | --------------- | -------- | --------- | ------------ |
| ONT-HG001                  | ONT R9.4.1            | HG001           | 55       | v3.3.2    | Guppy4.2.2   |
| ONT-HG001\_GP2.3.8         | ONT R9.4.1            | HG001           | 34       | v3.3.2    | Guppy2.3.8   |
| ONT-HG001\_GP2.3.8-4.2.2 | ONT R9.4.1            | HG001           | 45       | v3.3.2    | Guppy (2.3.8 + 4.2.2)   |
| ONT-HG001-4\_GP4.2.2     | ONT R9.4.1            | HG001-4         | 69       | v3.3.2 (HG001) + v4.2.1 (HG002-4)| Guppy4.2.2|
| ONT-HG002                  | ONT R9.4.1            | HG002           | 47       | v4.2.1    | Guppy4.2.2   |
| ONT-HG002\_GP4.2.2\_v3.3.2 | ONT R9.4.1            | HG002           | 47       | v3.3.2    | Guppy4.2.2   |
| ONT-HG002\_GP2.3.4\_v3.3.2 | ONT R9.4.1            | HG002           | 53       | v3.3.2    | Guppy2.3.4   |
| ONT-HG002\_GP2.3.4\_v4.2.1 | ONT R9.4.1            | HG002           | 53       | v4.2.1    | Guppy2.3.4   |
| ONT-HG002\_bonito          | ONT R9.4.1            | HG002 (chr1-21) | 51       | v4.2.1    | Bonito v0.30 |
| ONT-HG002\_r10.3           | ONT R10.3             | HG002 (chr1-21) | 32       | v4.2.1    | Guppy4.0.11  |
| CCS-HG001                  | PacBio CCS            | HG001           | 57       | v3.3.2    | \-           |
| CCS-HG002                  | PacBio CCS            | HG002           | 56       | v4.2.1    | \-           |
| CCS-HG001-4                | PacBio CCS            | HG001-4         | 55       | v3.3.2 (HG001) + v4.2.1 (HG002-4)| Guppy4.2.2    | \-           |
| CLR-HG002                  | PacBio CLR            | HG002           | 58       | v4.2.1    | \-           |
| NanoCaller1                | ONT R9.4.1            | HG001           | 34       | v3.3.2    | Guppy2.3.8   |
| NanoCaller2                | ONT R9.4.1            | HG002           | 53       | v3.3.2    | Guppy2.3.4   |
| NanoCaller3                | PacBio CLR            | HG003           | 28       | v3.3.2    | \-           |


#### Indel Models
| Model Name | Sequencing Technology | Genome | Coverage | Benchmark | Basecaller |
| ------------ | --------------------- | ------ | -------- | --------- | ---------- |
| ONT-HG001    | ONT R9.4.1            | HG001  | 55       | v3.3.2    | Guppy4.2.2 |
| ONT-HG002    | ONT R9.4.1            | HG002  | 47       | v4.2.1    | Guppy4.2.2 |
| CCS-HG001    | PacBio CCS            | HG001  | 57       | v3.3.2    | \-         |
| CCS-HG002    | PacBio CCS            | HG002  | 56       | v4.2.1    | \-         |
| NanoCaller1  | ONT R9.4.1            | HG001  | 34       | v3.3.2    | Guppy2.3.8 |
| NanoCaller3  | PacBio CCS            | HG001  | 29       | v3.3.2    | \-         |

## Citing NanoCaller
Please cite: Ahsan, M.U., Liu, Q., Fang, L. et al. NanoCaller for accurate detection of SNPs and indels in difficult-to-map regions from long-read sequencing by haplotype-aware deep neural networks. Genome Biol 22, 261 (2021). https://doi.org/10.1186/s13059-021-02472-2.
