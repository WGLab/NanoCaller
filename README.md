# NanoCaller
NanoCaller is a computational method that integrates long reads in deep convolutional neural network for the detection of SNPs/indels from long-read sequencing data. NanoCaller uses long-range haplotype structure to generate predictions for each SNP candidate variant site by considering pileup information of other candidate sites sharing reads. Subsequently, it performs read phasing, and carries out local realignment of each set of phased reads and the set of all reads for each indel candidate variant site to generate indel calling, and then creates consensus sequences for indel sequence prediction.

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
