# NanoCaller
NanoCaller is a computational method that integrates long reads in deep convolutional neural network for the detection of SNPs/indels from long-read sequencing data. NanoCaller uses long-range haplotype structure to generate predictions for each SNP candidate variant site by considering pileup information of other candidate sites sharing reads. Subsequently, it performs read phasing, and carries out local realignment of each set of phased reads and the set of all reads for each indel candidate variant site to generate indel calling, and then creates consensus sequences for indel sequence prediction.

## Latest Updates
_**v0.4.0** (June 2 2021)_ : Added NanoCaller models trained on ONT reads basecalled with Guppy v4.2.2 and Bonito v0.30, as well as R10.3 reads. Added new NanoCaller models trained with long CCS reads (15-20kb library selection). Improved indel calling with rolling window for candidate selection which helps with indels in low complexity regions.

## Citing NanoCaller
Please cite: Ahsan, Umair and Liu, Qian and Wang, Kai. NanoCaller for accurate detection of SNPs and small indels from long-read sequencing by deep neural networks. bioRxiv 2019.12.29.890418; doi: https://doi.org/10.1101/2019.12.29.890418

## Installation
NanoCaller can be installed using Docker or Conda. Please refer to [Installation](docs/Install.md) for instructions regarding installing NanoCaller.

## Usage
General usage of NanoCaller is described in [Usage](docs/Usage.md). For a comprehensive case study of variant calling on Nanopore reads, see [ONT Case Study](docs/ONT%20Case%20Study.md), where we describe end-to-end variant calling pipeline for using NanoCaller, where we start with aligning FASTQ files of HG002, calls variants using NanoCaller, and evaluate performances on various genomic regions.

## Example
An example of NanoCaller usage is provided in [sample](sample). The results are stored in [test output](sample/test_run) and were created using the following command:

`python ../scripts/NanoCaller.py -bam HG002.nanopore.chr22.sample.bam -p ont -o test_run -chrom chr22 -start 20000000 -end 21000000 -ref chr22_ref.fa -cpu 4 > log`

which is also in the file [sample_call](sample/sample_call). This example should take about 10-15 minutes to run.

## Well-trained models
The well-trained models for [ONT](https://github.com/WGLab/NanoCaller/tree/master/scripts/release_data/ONT_models) data, [CLR](https://github.com/WGLab/NanoCaller/tree/master/scripts/release_data/clr_models) data and [HIFI](https://github.com/WGLab/NanoCaller/tree/master/scripts/release_data/hifi_models) data can be found [here](https://github.com/WGLab/NanoCaller/tree/master/scripts/release_data). These models are trained on chr1-22 of the genomes they are trained on, unless mentioned othewise.

You can specify SNP and indel models using `--snp_model` and `--indel_model` parameters with a model name from tables below. For instance, if you want to use 'ONT-HG002\_bonito' SNP model and 'ONT-HG002' indel model, use the following command:

`python NanoCaller.py --snp_model ONT-HG002_bonito--indel_model ONT-HG002`



Feel free to share your models if you use NanoCallers to train other models. 

#### SNP Models

| Model Name                 | Sequencing Technology | Genome          | Coverage | Benchmark | Basecaller   |
| -------------------------- | --------------------- | --------------- | -------- | --------- | ------------ |
| ONT-HG001                  | ONT R9.4.1            | HG001           | 55       | v3.3.2    | Guppy4.2.2   |
| ONT-HG001\_GP2.3.8         | ONT R9.4.1            | HG001           | 34       | v3.3.2    | Guppy2.3.8   |
| ONT-HG001_guppy2.3.8_guppy4.2.2 | ONT R9.4.1            | HG001           | 45       | v3.3.2    | Guppy (2.3.8 + 4.2.2)   |
| ONT-HG002                  | ONT R9.4.1            | HG002           | 47       | v4.2.1    | Guppy4.2.2   |
| ONT-HG002\_GP4.2.2\_v3.3.2 | ONT R9.4.1            | HG002           | 47       | v3.3.2    | Guppy4.2.2   |
| ONT-HG002\_GP2.3.4\_v3.3.2 | ONT R9.4.1            | HG002           | 53       | v3.3.2    | Guppy2.3.4   |
| ONT-HG002_guppy2.3.4_giab-4.2.1 | ONT R9.4.1            | HG002           | 53       | v4.2.1    | Guppy2.3.4   |
| ONT-HG002\_bonito          | ONT R9.4.1            | HG002 (chr1-21) | 51       | v4.2.1    | Bonito v0.30 |
| ONT-HG002\_r10.3           | ONT R10.3             | HG002 (chr1-21) | 32       | v4.2.1    | Guppy4.0.11  |
| ONT-HG001-4_guppy4.2.2     | ONT R9.4.1            | HG001-4         | 69       | v3.3.2 (HG001) + v4.2.1 (HG002-4)| Guppy4.2.2|
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
