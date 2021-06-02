# Usage

### Using Docker to run NanoCaller

Please check the [NanoCaller Docker Hub repository](https://hub.docker.com/repository/docker/genomicslab/nanocaller) for the most up to date version of NanoCaller docker image.

#### Whole Genome Variant Calling

For whole genome variant calling, or calling variants on several chromosomes, use `NanoCaller_WGS.py` to call variants. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. 
```
VERSION=0.4.0
docker run -itd \
-v 'YOUR_INPUT_DIR':'/input/' \
-v 'YOUR_WORKING_DIR':'/output/' \
genomicslab/nanocaller:${VERSION} \
python NanoCaller_WGS.py \
-bam /input/YOUR_BAM \
-ref /input/YOUR_REF \
-o /output \
-cpu NUMBER_OF_CPUS_TO_USE
-p PRESET
```

#### Single Chromosome Variant Calling

For calling variants on single chromosomes, use `NanoCaller.py` to call variants. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. 
```
VERSION=0.4.0
docker run -itd \
-v 'YOUR_INPUT_DIR':'/input/' \
-v 'YOUR_WORKING_DIR':'/output/' \
genomicslab/nanocaller:${VERSION} \
python NanoCaller.py \
-chrom CHROMOSOME \
-bam /input/YOUR_BAM \
-ref /input/YOUR_REF \
-o /output \
-cpu NUMBER_OF_CPUS_TO_USE 
-p PRESET
```



### Running NanoCaller without Docker
#### Whole Genome Variant Calling

For whole genome variant calling, or calling variants on several chromosomes, use `NanoCaller_WGS.py` to call variants. Assuming NanoCaller respository is cloned and located at `PATH_TO_NANOCALLER_REPOSITORY` directory, run the following command.
```
python PATH_TO_NANOCALLER_REPOSITORY/scripts/NanoCaller_WGS.py 
-bam YOUR_BAM \
-ref YOUR_REF \
-o OUTPUT_DIRECTORY \
-cpu NUMBER_OF_CPUS_TO_USE
-p PRESET
```

Another way to run NanoCaller for whole genome variant calling is to use `NanoCaller.py` or `NanoCaller_WGS.py` on each chromosome separately by setting `chrom` argument. This option is suitable if you have a large computing cluster with a lot of computing resources. For instance, on a Sun Grid Engine, you can submit a separate job for each chromosome like this, using 16 CPUs per job:

```
for i in {1..22};do echo "python PATH_TO_NANOCALLER_REPOSITORY/scripts/NanoCaller.py
-chrom chr$i
-bam YOUR_BAM \
-ref YOUR_REF \
-o OUTPUT_DIRECTORY \
-cpu 16
-p PRESET" |qsub -V -cwd -pe smp 16 -N chr$i -e chr$i.log -o chr$i.log; done
```



#### Single Chromosome Variant Calling

For calling variants on single chromosomes, use `NanoCaller.py` to call variants.

 Assuming NanoCaller respository is cloned and located at `PATH_TO_NANOCALLER_REPOSITORY` directory, run the following command.
```
python PATH_TO_NANOCALLER_REPOSITORY/scripts/NanoCaller_WGS.py 
-chrom CHROMOSOME \
-bam YOUR_BAM \
-ref YOUR_REF \
-o OUTPUT_DIRECTORY \
-cpu NUMBER_OF_CPUS_TO_USE 
-p PRESET
```

## General Usage Options

### Single Chromosome Variant Calling

### Whole Genome Variant Calling
```
usage: NanoCaller_WGS.py [-h] [-mode MODE] [-seq SEQUENCING] [-cpu CPU]
                         [-mincov MINCOV] [-maxcov MAXCOV] [-keep_bam]
                         [-o OUTPUT] [-prefix PREFIX] [-sample SAMPLE]
                         [-chrom [CHROM [CHROM ...]]]
                         [-include_bed INCLUDE_BED] [-exclude_bed EXCLUDE_BED]
                         [-wgs_contigs_type WGS_CONTIGS_TYPE] [-p PRESET] -bam
                         BAM -ref REF [-snp_model SNP_MODEL]
                         [-min_allele_freq MIN_ALLELE_FREQ]
                         [-min_nbr_sites MIN_NBR_SITES]
                         [-nbr_t NEIGHBOR_THRESHOLD] [-sup]
                         [-indel_model INDEL_MODEL] [-ins_t INS_THRESHOLD]
                         [-del_t DEL_THRESHOLD] [-win_size WIN_SIZE]
                         [-small_win_size SMALL_WIN_SIZE] [-phase_bam]
                         [-enable_whatshap]

optional arguments:
  -h, --help            show this help message and exit

Required Arguments:
  -bam BAM, --bam BAM   Bam file, should be phased if 'indel' mode is selected
                        (default: None)
  -ref REF, --ref REF   Reference genome file with .fai index (default: None)

Preset:
  -p PRESET, --preset PRESET
                        Apply recommended preset values for SNP and Indel
                        calling parameters, options are 'ont', 'ul_ont',
                        'ul_ont_extreme', 'ccs' and 'clr' (default: None)

Configurations:
  -mode MODE, --mode MODE
                        NanoCaller mode to run, options are 'snps',
                        'snps_unphased', 'indels' and 'both'. 'snps_unphased'
                        mode quits NanoCaller without using WhatsHap for
                        phasing. (default: both)
  -seq SEQUENCING, --sequencing SEQUENCING
                        Sequencing type, options are 'ont', 'ul_ont',
                        'ul_ont_extreme', and 'pacbio' (default: ont)
  -cpu CPU, --cpu CPU   Number of CPUs to use (default: 1)
  -mincov MINCOV, --mincov MINCOV
                        Minimum coverage to call a variant (default: 8)
  -maxcov MAXCOV, --maxcov MAXCOV
                        Maximum coverage of reads to use. If sequencing depth
                        at a candidate site exceeds maxcov then reads are
                        downsampled. (default: 160)

Variant Calling Regions:
  -chrom [CHROM [CHROM ...]], --chrom [CHROM [CHROM ...]]
                        A space/whitespace separated list of contigs, e.g.
                        chr3 chr6 chr22. (default: None)
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
  -wgs_contigs_type WGS_CONTIGS_TYPE, --wgs_contigs_type WGS_CONTIGS_TYPE
                        Options are "with_chr", "without_chr" and "all",
                        "with_chr" option will assume human genome and run
                        NanoCaller on chr1-22, "without_chr" will run on
                        chromosomes 1-22 if the BAM and reference genome files
                        use chromosome names without "chr". "all" option will
                        run NanoCaller on each contig present in reference
                        genome FASTA file. (default: with_chr)

SNP Calling:
  -snp_model SNP_MODEL, --snp_model SNP_MODEL
                        NanoCaller SNP model to be used (default: ONT-HG002)
  -min_allele_freq MIN_ALLELE_FREQ, --min_allele_freq MIN_ALLELE_FREQ
                        minimum alternative allele frequency (default: 0.15)
  -min_nbr_sites MIN_NBR_SITES, --min_nbr_sites MIN_NBR_SITES
                        minimum number of nbr sites (default: 1)
  -nbr_t NEIGHBOR_THRESHOLD, --neighbor_threshold NEIGHBOR_THRESHOLD
                        SNP neighboring site thresholds with lower and upper
                        bounds seperated by comma, for Nanopore reads
                        '0.4,0.6' is recommended, for PacBio CCS anc CLR reads
                        '0.3,0.7' and '0.3,0.6' are recommended respectively
                        (default: 0.4,0.6)
  -sup, --supplementary
                        Use supplementary reads (default: False)

Indel Calling:
  -indel_model INDEL_MODEL, --indel_model INDEL_MODEL
                        NanoCaller indel model to be used (default: ONT-HG002)
  -ins_t INS_THRESHOLD, --ins_threshold INS_THRESHOLD
                        Insertion Threshold (default: 0.4)
  -del_t DEL_THRESHOLD, --del_threshold DEL_THRESHOLD
                        Deletion Threshold (default: 0.6)
  -win_size WIN_SIZE, --win_size WIN_SIZE
                        Size of the sliding window in which the number of
                        indels is counted to determine indel candidate site.
                        Only indels longer than 2bp are counted in this
                        window. Larger window size can increase recall, but
                        use a maximum of 50 only (default: 40)
  -small_win_size SMALL_WIN_SIZE, --small_win_size SMALL_WIN_SIZE
                        Size of the sliding window in which indel frequency is
                        determined for small indels (default: 4)

Output Options:
  -keep_bam, --keep_bam
                        Keep phased bam files. (default: False)
  -o OUTPUT, --output OUTPUT
                        VCF output path, default is current working directory
                        (default: None)
  -prefix PREFIX, --prefix PREFIX
                        VCF file prefix (default: variant_calls)
  -sample SAMPLE, --sample SAMPLE
                        VCF file sample name (default: SAMPLE)

Phasing:
  -phase_bam, --phase_bam
                        Phase bam files if snps mode is selected. This will
                        phase bam file without indel calling. (default: False)
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

```
## Understanding NanoCaller Output
Depending upon which mode is run, NanoCaller will produce the following files:

- PREFIX.snps.vcf.gz contains unphased SNP calls made by NanoCaller using a deep learning model. NanoCaller modes that produce this file are: `snps_unphased`, `snps` and `both`.
- PREFIX.snps.phased.vcf.gz contains SNP calls from PREFIX.snps.vcf.gz that are phase with WhatsHap. By default they have the same genotype as in the PREFIX.snps.vcf.gz file, unless `--enable_whatshap` flag is set which can allow WhatsHap to change genotypes. NanoCaller modes that produce this file are: `snps` and `both`.
- PREFIX.indels.vcf.gz contains indel calls made by NanoCaller using multiple sequence alignment. Some of these calls might be indels combined with nearby substitutions or multi-nucleotide substitutions. NanoCaller modes that produce this file are: `indels` and `both`.
- PREFIX.final.vcf.gz contains SNP calls from PREFIX.snps.phased.vcf.gz and indel calls from PREFIX.indels.vcf.gz. NanoCaller mode that produce this file is: `both`.

## Preset Options
Users can select recommended NanoCaller settings for various sequencing types by using preset option `-p` or `--preset`. There are five presets, which are described below with their equivalent parameter settings. Any parameters included in the presets can be overwritten by specifying the parameter.

### ONT reads
Options are: `ont`, `ul_ont` and `ul_ont_extreme`. These presets use ONT trained models for SNP and indell calling, and only differ in `--sequencing` argument. With `ont`, `ul_ont` and `ul_ont_extreme` options, NanoCaller uses haplotype information from up to 50kbp, 100kbp and 300kbp bases away, respectively, for SNP calling. There is no difference in indel calling between these options. Therefore, use `ul_ont` if you have a significant number of reads that are up to 100kbp long, and `ul_ont_extreme` if you have a significant number of reads that are up to 300kbp long.

Preset `ont` is equivalent to:
```
--sequencing ont
--snp_model ONT-HG002
--indel_model ONT-HG002
--neighbor_threshold '0.4,0.6'
--ins_threshold 0.4
--del_threshold 0.6
```

Preset `ul_ont` is equivalent to:
```
--sequencing ul_ont
--snp_model ONT-HG002
--indel_model ONT-HG002
--neighbor_threshold '0.4,0.6'
--ins_threshold 0.4
--del_threshold 0.6
```

Preset `ul_ont_extreme` is equivalent to:
```
--sequencing ul_ont_extreme
--snp_model ONT-HG002
--indel_model ONT-HG002
--neighbor_threshold '0.4,0.6'
--ins_threshold 0.4
--del_threshold 0.6
```

### PacBio reads:
Options are `ccs` and `clr`, for PacBio CCS (HiFi) and CLR reads respectively.

Preset `ccs` is equivalent to:
```
--sequencing pacbio
--snp_model CCS-HG002
--indel_model CCS-HG002
--neighbor_threshold '0.3,0.7'
--ins_threshold 0.4
--del_threshold
--enable_whatshap
```

Preset `clr` is equivalent to:
```
--sequencing pacbio
--snp_model CLR-HG002
--indel_model ONT-HG002
--neighbor_threshold '0.3,0.6'
--ins_threshold 0.6
--del_threshold 0.6
--win_size 10
--small_win_size 2
--enable_whatshap
```


## Important Considerations

Some important options to keep in mind when using NanoCaller:
- `--sequencing` argument is important to set because NanoCaller has slightly different settings for generating inputs for each type of sequencing. If you use a preset option `--preset`, then `--sequencing`  is set automatically.
- `--exclude_bed` argument can be used to speed up NanoCaller. For instance, by setting it to `hg38`, you can tell NanoCaller to skip centromere and telomere regions which have incredibly high number of candidate variants due to poor alignment.
- `--mode` argument can be used to select which types of variants you want to call. The fastest option is `snps_unphased` which only makes SNP calls.
- `--phase_bam` flag can be set to phase reads if you are only interested in calling SNPs, and phasing input BAM file with WhatsHap. This argument only works when you use `--mode snps`, and it is recommended to use it with NanoCaller's single variant calling routine using `NanoCaller.py` to phase reads from the entire chromosome together. If you use `--mode both` then reads are automatically phased for indel calling and you can select to keep phased BAM files them using `--keep_bam`. 
- `-nbr_t` option is sensitive to sequencing type so choose this accordingly.
- `-ins_t` and `del_t` are insertion and deletion frequency thresholds are per haplotype. Default values are slightly higher due to high error in ONT reads, but these thresholds can be lowered for CCS reads.
- CLR reads have incredibly low insertion and deletion freqencies in a pileup due highly variable placement of indels by aligners. To detect indels on CLR reads, you would need to set low frequency threhsolds, which leads to a huge increase in runtime.
- `--keep_bam` flag can be used to save phased BAM files created by NanoCaller for indel calling. By default, we delete these BAM files in order to not use up too much storage.
- `--enable_whatshap` flag can improve SNP calling performance for PacBio reads.

We recommend using `NanoCaller_WGS.py` for whole genome variant calling, and `NanoCaller.py` for single chromosome variant calling, although `NanoCaller_WGS.py` can also be used for single chromosome. `NanoCaller_WGS.py` breaks genome into 10Mb chunks, uses GNU parallel to run `NanoCaller.py` on each chunk independently using 1 CPU, and then combines the results at the end. We do not break genome into chunks smaller than 10Mb so that phasing can be done accurately. Lets say you can use 20 CPUs, and if you use `NanoCaller_WGS.py` for whole human genome of 3000Mb, NanoCaller will create 300 jobs for 20 CPUs to run. But if you use it for a chromosome of size 100Mb, then only 10 jobs will be created for 20 CPUs, and using more than 10 CPUs is not going to improve runtime and the rest of the CPUs will be left idle. If you run `NanoCaller.py` on the same 100Mb chromosome with 20 CPUs, `NanoCaller.py` will use python's multiprocessing module to generate features for SNPs in chunks of 200Kb using 20 CPUs by running 500 jobs, use only a single CPU to run WhatsHap for phasing the entire chromosome, followed by generating features for indels in chunks of 50Kb using 20 CPUs by running 2000 jobs. This method is more aggressive in terms of utilizing computing resources by breaking the chromosome into very small chunks for feature generation, creates fewer intermediate files and folder, but can have bottleneck issues with using WhatsHap with 1 CPU only.
