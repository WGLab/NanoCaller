# Usage

### Using Docker to run NanoCaller
#### Whole Genome Variant Calling

For whole genome variant calling, or calling variants on several chromosomes, use `NanoCaller_WGS.py` to call variants. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. 
```
VERSION=0.3.0
docker run -itd \
-v 'YOUR_INPUT_DIR':'/input/' \
-v 'YOUR_WORKING_DIR':'/output/' \
genomicslab/nanocaller:${VERSION} \
python NanoCaller_WGS.py \
-bam /input/YOUR_BAM \
-ref /input/YOUR_REF \
-prefix PREFIX \
-o /output \
-cpu NUMBER_OF_CPUS_TO_USE
-seq SEQUENCING_TYPE
```

#### Single Chromosome Variant Calling

For calling variants on single chromosomes, use `NanoCaller.py` to call variants. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. 
```
VERSION=0.3.0
docker run -itd \
-v 'YOUR_INPUT_DIR':'/input/' \
-v 'YOUR_WORKING_DIR':'/output/' \
genomicslab/nanocaller:${VERSION} \
python NanoCaller.py \
-chrom CHROMOSOME \
-bam /input/YOUR_BAM \
-ref /input/YOUR_REF \
-prefix PREFIX \
-o /output \
-cpu NUMBER_OF_CPUS_TO_USE 
-seq SEQUENCING_TYPE
```



### Running NanoCaller without Docker
#### Whole Genome Variant Calling

For whole genome variant calling, or calling variants on several chromosomes, use `NanoCaller_WGS.py` to call variants. Assuming NanoCaller respository is cloned and located at `PATH_TO_NANOCALLER_REPOSITORY` directory, run the following command.
```
python PATH_TO_NANOCALLER_REPOSITORY/scripts/NanoCaller_WGS.py 
-bam YOUR_BAM \
-ref YOUR_REF \
-prefix PREFIX \
-o OUTPUT_DIRECTORY \
-cpu NUMBER_OF_CPUS_TO_USE
-seq SEQUENCING_TYPE
```

Another way to run NanoCaller for whole genome variant calling is to use `NanoCaller.py` or `NanoCaller_WGS.py` on each chromosome separately by setting `chrom` argument. This option is suitable if you have a large computing cluster with a lot of computing resources. For instance, on a Sun Grid Engine, you can submit a separate job for each chromosome like this, using 16 CPUs per job:

```
for i in {1..22};do echo "python PATH_TO_NANOCALLER_REPOSITORY/scripts/NanoCaller.py
-chrom chr$i
-bam YOUR_BAM \
-ref YOUR_REF \
-prefix PREFIX \
-o OUTPUT_DIRECTORY \
-cpu 16
-seq SEQUENCING_TYPE" |qsub -V -cwd -pe smp 16 -N chr$i -e chr$i.log -o chr$i.log; done
```



#### Single Chromosome Variant Calling

For calling variants on single chromosomes, use `NanoCaller.py` to call variants.

 Assuming NanoCaller respository is cloned and located at `PATH_TO_NANOCALLER_REPOSITORY` directory, run the following command.
```
python PATH_TO_NANOCALLER_REPOSITORY/scripts/NanoCaller_WGS.py 
-chrom CHROMOSOME \
-bam YOUR_BAM \
-ref YOUR_REF \
-prefix PREFIX \
-o OUTPUT_DIRECTORY \
-cpu NUMBER_OF_CPUS_TO_USE 
-seq SEQUENCING_TYPE
```

## General Usage Options
```
usage: NanoCaller_WGS.py [-h] [-mode MODE] [-seq SEQUENCING] [-model MODEL]
                         [-o OUTPUT] [-chrom CHROM] [-cpu CPU]
                         [-min_allele_freq MIN_ALLELE_FREQ]
                         [-min_nbr_sites MIN_NBR_SITES] -bam BAM -ref REF
                         -prefix PREFIX [-include_bed INCLUDE_BED]
                         [-exclude_bed EXCLUDE_BED] [-sample SAMPLE] [-sup]
                         [-mincov MINCOV] [-maxcov MAXCOV] [-start START]
                         [-end END] [-nbr_t NEIGHBOR_THRESHOLD]
                         [-ins_t INS_THRESHOLD] [-del_t DEL_THRESHOLD]
                         [-enable_whatshap] [-keep_bam]
                         [-wgs_contigs_type WGS_CONTIGS_TYPE]

optional arguments:
  -h, --help            show this help message and exit
  -mode MODE, --mode MODE
                        Testing mode, options are 'snps', 'snps_unphased',
                        'indels' and 'both'. 'snps_unphased' mode quits
                        NanoCaller without using WhatsHap for phasing.
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
  -o OUTPUT, --output OUTPUT
                        VCF output path, default is current working directory
                        (default: None)
  -chrom CHROM, --chrom CHROM
                        A space/whitespace separated list of contigs in
                        quotation marks e.g. "chr3 chr6 chr22" . (default:
                        None)
  -cpu CPU, --cpu CPU   CPUs (default: 1)
  -min_allele_freq MIN_ALLELE_FREQ, --min_allele_freq MIN_ALLELE_FREQ
                        minimum alternative allele frequency (default: 0.15)
  -min_nbr_sites MIN_NBR_SITES, --min_nbr_sites MIN_NBR_SITES
                        minimum number of nbr sites (default: 1)
  -include_bed INCLUDE_BED, --include_bed INCLUDE_BED
                        Only call variants inside the intervals specified in
                        the bgzipped and tabix indexed BED file. If any other
                        flags are used to specify a region, intersect the
                        region with intervals in the BED file, e.g. if -chom
                        chr1 -start 10000000 -end 20000000 flags are set, call
                        variants inside the intervals specified by the BED
                        file that overlap with chr1:10000000-20000000. Same
                        goes for the case when whole genome variant calling
                        flag is set. (default: None)
  -exclude_bed EXCLUDE_BED, --exclude_bed EXCLUDE_BED
                        Path to bgzipped and tabix indexed BED file containing
                        intervals to ignore for variant calling. BED files of
                        centromere and telomere regions for the following
                        genomes are included in NanoCaller: hg38, hg19, mm10
                        and mm39. To use these BED files use one of the
                        following options: ['hg38', 'hg19', 'mm10', 'mm39'].
                        (default: None)
  -sample SAMPLE, --sample SAMPLE
                        VCF file sample name (default: SAMPLE)
  -sup, --supplementary
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
                        '0.4,0.6' is recommended, for PacBio CCS anc CLR reads
                        '0.3,0.7' and '0.3,0.6' are recommended respectively
                        (default: 0.4,0.6)

  -ins_t INS_THRESHOLD, --ins_threshold INS_THRESHOLD
                        Insertion Threshold (default: 0.4)
  -del_t DEL_THRESHOLD, --del_threshold DEL_THRESHOLD
                        Deletion Threshold (default: 0.6)
  -enable_whatshap, --enable_whatshap
                        Allow WhatsHap to change SNP genotypes when phasing
                        using --distrust-genotypes and --include-homozygous
                        flags (this is not the same as regenotyping),
                        considerably increasing the time needed for phasing.
                        It has a negligible effect on SNP calling accuracy for
                        Nanopore reads, but may make a small improvement for
                        PacBio reads. By default WhatsHap will only phase SNP
                        calls produced by NanoCaller, but not change their
                        genotypes. (default: False)
  -keep_bam, --keep_bam
                        Keep phased bam files. (default: False)
  -wgs_contigs_type WGS_CONTIGS_TYPE, --wgs_contigs_type WGS_CONTIGS_TYPE
                        Options are "with_chr", "without_chr" and "all",
                        "with_chr" option will assume human genome and run
                        NanoCaller on chr1-22, "without_chr" will run on
                        chromosomes 1-22 if the BAM and reference genome files
                        use chromosome names without "chr". "all" option will
                        run NanoCaller on each contig present in reference
                        genome FASTA file. (default: with_chr)

Required arguments:
  -bam BAM, --bam BAM   Bam file, should be phased if 'indel' mode is selected
                        (default: None)
  -ref REF, --ref REF   reference genome file with .fai index (default: None)
  -prefix PREFIX, --prefix PREFIX
                        VCF file prefix (default: None)
```

## Usage Cases and Suggestions

Some important options to keep in mind when using NanoCaller:
- `seq` argument is important to set to either `ont` or `pacbio` because NanoCaller has slightly different settings for generating inputs for each type of sequencing.
- `--exclude_bed` argument can be used to speed up NanoCaller. For instance, by setting it to `hg38`, you can tell NanoCaller to skip centromere and telomere regions which have incredibly high number of candidate variants due to poor alignment.
- `model` argument can be used to choose which SNP calling model you want to use. By default we use `NanoCaller1` which is trained on HG001 ONT reads, and it works very well across genomes and sequencing technologies. We do not recommend using `NanoCaller3` model trained on PacBio data to call variants on ONT reads.
- `mode` argument can be used to select which types of variants you want to call. The fastest option is `snps_unphased` which only make SNP calls.
- `nbr_t` option is sensitive to sequencing type so choose this accordingly.
- `ins_t` and `del_t` are insertion and deletion frequency thresholds are per haplotype. Default values are slightly higher due to high error in ONT reads, but these thresholds can be lowered for CCS reads.
- We do not recommend using NanoCaller to call indels on PacBio CLR reads. CLR reads have incredibly low insertion and deletion freqencies in a pileup due highly variable placement of indels by aligners. To detect indels on CLR reads, you would need to set low frequency threhsolds, which leads to a huge increase in runtime.
- `keep_bam` flag can be used to save phased BAM files created by NanoCaller for indel calling. By default, we delete these BAM files in order to not use up too much storage.

We recommend using `NanoCaller_WGS.py` for whole genome variant calling, and `NanoCaller.py` for single chromosome variant calling, although `NanoCaller_WGS.py` can also be used for single chromosome. `NanoCaller_WGS.py` breaks genome into 10Mb chunks, uses GNU parallel to run `NanoCaller.py` on each chunk independently using 1 CPU, and then combines the results at the end. We do not break genome into chunks smaller than 10Mb so that phasing can be done accurately. Lets say you can use 20 CPUs, and if you use `NanoCaller_WGS.py` for whole human genome of 3000Mb, NanoCaller will create 300 jobs for 20 CPUs to run. But if you use it for a chromosome of size 100Mb, then only 10 jobs will be created for 20 CPUs, and using more than 10 CPUs is not going to improve runtime and the rest of the CPUs will be left idle. If you run `NanoCaller.py` on the same 100Mb chromosome with 20 CPUs, `NanoCaller.py` will use python's multiprocessing module to generate features for SNPs in chunks of 200Kb using 20 CPUs by running 500 jobs, use only a single CPU to run WhatsHap for phasing the entire chromosome, followed by generating features for indels in chunks of 50Kb using 20 CPUs by running 2000 jobs. This method is more aggressive in terms of utilizing computing resources by breaking the chromosome into very small chunks for feature generation, creates fewer intermediate files and folder, but can have bottleneck issues with using WhatsHap with 1 CPU only.
