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
After Anaconda is installed successfully, add the following channels using the commands:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

After adding these channels, simply run the following commands to install all NanoCaller dependencies into a conda environment 'nanocaller_env':

```
conda create -n nanocaller_env -c bioconda bcftools biopython  muscle numpy pysam python=3.6.8 rtg-tools samtools=1.10 tensorflow=1.13.* whatshap=1.0 vcflib

conda activate nanocaller_env
```
Clone the NanoCaller github repository
```
git clone https://github.com/WGLab/NanoCaller.git
cd NanoCaller
```

## Usage
```
python PATH_TO_NANOCALLER_REPOSITORY/scripts/NanoCaller.py [-h] [-mode MODE] 
		     [-seq SEQUENCING] [-model MODEL]
                     [-vcf VCF] -chrom CHROM [-cpu CPU]
                     [-min_allele_freq MIN_ALLELE_FREQ]
                     [-min_nbr_sites MIN_NBR_SITES] -bam BAM -ref REF -prefix
                     PREFIX [-sample SAMPLE] [-sup SUPPLEMENTARY]
                     [-mincov MINCOV] [-maxcov MAXCOV] [-start START]
                     [-end END] [-nbr_t NEIGHBOR_THRESHOLD]
                     [-ins_t INS_THRESHOLD] [-del_t DEL_THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -mode MODE, --mode MODE
                        Testing mode, options are 'snps', 'indels' and 'both'
                        (default: both)
  -seq SEQUENCING, --sequencing SEQUENCING
                        Sequencing type, options are 'ont' and 'pacbio'
                        (default: ont)
  -model MODEL, --model MODEL
                        NanoCaller SNP model to be used, options are
                        'NanoCaller1' (trained on HG001 Nanopore reads),
                        'NanoCaller2' (trained on HG002 Nanopore reads) and
                        'NanoCaller3' (trained on HG003 PacBio reads)
                        (default: NanoCaller1)
  -vcf VCF, --vcf VCF   VCF output path, default is current working directory
                        (default: None)
  -cpu CPU, --cpu CPU   CPUs (default: 1)
  -min_allele_freq MIN_ALLELE_FREQ, --min_allele_freq MIN_ALLELE_FREQ
                        minimum alternative allele frequency (default: 0.15)
  -min_nbr_sites MIN_NBR_SITES, --min_nbr_sites MIN_NBR_SITES
                        minimum number of nbr sites (default: 1)
  -sample SAMPLE, --sample SAMPLE
                        VCF file sample name (default: SAMPLE)
  -sup SUPPLEMENTARY, --supplementary SUPPLEMENTARY
                        Use supplementary reads (default: False)
  -mincov MINCOV, --mincov MINCOV
                        min coverage (default: 8)
  -maxcov MAXCOV, --maxcov MAXCOV
                        max coverage (default: 160)
  -start START, --start START
                        start, default is 1 (default: None)
  -end END, --end END   end, default is the end of contig (default: None)
  -nbr_t NEIGHBOR_THRESHOLD, --neighbor_threshold NEIGHBOR_THRESHOLD
                        SNP neighboring site thresholds with lower and upper
                        bounds seperated by comma, for Nanopore reads
                        '0.4,0.6' is recommended and for PacBio reads
                        '0.3,0.7' is recommended (default: 0.4,0.6)
  -ins_t INS_THRESHOLD, --ins_threshold INS_THRESHOLD
                        Insertion Threshold (default: 0.4)
  -del_t DEL_THRESHOLD, --del_threshold DEL_THRESHOLD
                        Deletion Threshold (default: 0.6)

Required arguments:
  -chrom CHROM, --chrom CHROM
                        Chromosome (default: None)
  -bam BAM, --bam BAM   Bam file, should be phased if 'indel' mode is selected
                        (default: None)
  -ref REF, --ref REF   reference genome file with .fai index (default: None)
  -prefix PREFIX, --prefix PREFIX
                        VCF file prefix (default: None)
```
## Example
An example of NanoCaller usage is provided in [sample](sample). The results are stored in [test output](sample/test_run) and were created using the following command:

`python ../scripts/NanoCaller.py -bam HG002.nanopore.chr22.sample.bam -mode both -seq ont -model NanoCaller1 -vcf test_run -chrom chr22 -start 20000000 -end 21000000 -ref chr22_ref.fa -prefix HG002.chr22.sample -cpu 1 > log`

which is also in the file [sample_call](sample/sample_call).
