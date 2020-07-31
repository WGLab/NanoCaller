# NanoCaller
NanoCaller is a computational method that integrates long reads in deep convolutional neural network for the detection of SNPs/indels from long-read sequencing data. NanoCaller uses long-range haplotype structure to generate predictions for each SNP candidate variant site by considering pileup information of other candidate sites sharing reads. Subsequently, it performs read phasing, and carries out local realignment of each set of phased reads and the set of all reads for each indel candidate variant site to generate indel calling, and then creates consensus sequences for indel sequence prediction.

## Citing NanoCaller
Please cite: Ahsan, Umair and Liu, Qian and Wang, Kai. NanoCaller for accurate detection of SNPs and small indels from long-read sequencing by deep neural networks. bioRxiv 2019.12.29.890418; doi: https://doi.org/10.1101/2019.12.29.890418

## Installation
First, install Miniconda, a minimal installation of Anaconda, which is much smaller and has a faster installation.
Note that this version is meant for Linux below, macOS and Windows have a different script:

```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Go through all the prompts (installation in `$HOME` is recommended).
After Anaconda is installed successfully, simply run:

```
conda create -n nanocaller_env -c bioconda bcftools=1.10.2 biopython=1.76 htslib=1.10.2  muscle=3.8.1551 numpy=1.18.5 pysam=0.15.3 python=3.6.8 rtg-tools=3.11 samtools=1.10 tensorflow=1.13.1 whatshap=1.0

git clone https://github.com/WGLab/NanoCaller.git
cd NanoCaller
conda activate nanocaller_env
```

## Usage
```
python NanoCaller.py -mode both       # what type of variants to call, option 'snps', 'indels and 'both 
-chrom chr1 \                         # contig name
-ref ref.fa \                         # FASTA reference file with .fai index
-bam alignments.bam \                 # BAM file with index (must be phased if 'indels' mode is chosen
-model NanoCaller1 \                  # model used for SNP calling
-bed regions.bed \                    # BED files 
-vcf sample_calls/ \                  # output directory
-cpu 4 \                              # number of CPUs to use
```
## Example
An example of NanoCaller usage is provided in [sample](sample). The results are stored in [test output](sample/test_run) and were created using the following command:

`python ../scripts/NanoCaller.py -bam HG002.nanopore.chr22.sample.bam -mode both -seq ont -model NanoCaller1 -vcf test_run -chrom chr22 -start 20000000 -end 21000000 -ref chr22_ref.fa -prefix HG002.chr22.sample -cpu 1 > log`

which is also in the file [sample_call](sample/sample_call).
