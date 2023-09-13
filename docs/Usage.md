# Usage

## NanoCaller Modes

NanoCaller can be run in three modes (`snps indels all`) that are specified using `--mode` parameter. 
- `snps`: this is the fastest run option which only procudes SNP calls abd can additionally output phased VCF and BAM files using WhatsHap if `--phase` parameter is also used
- `indels`: in this mode, NanoCaller expects a phased BAM file with HP tags to label read haplotypes. This mode will produce only indels shorter than 50bp, and may also produce any SNP calls that are within complex indel blocks consiting of several small variants in a short distance.
- `all`: the default mode which runs `snps` mode, phases BAM file and then runs `indels` mode on phased BAM files.

## Running NanoCaller with Docker

Please check the [NanoCaller Docker Hub repository](https://hub.docker.com/repository/docker/genomicslab/nanocaller) for the most up to date version of NanoCaller docker image.

NanoCaller can call variants from whole genome or several chromosomes using a single command `NanoCaller`. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. You can use `--regions` or `--bed` to specify regions for variant calling, otherwise all contigs in the input BAM file will be used. 

You can pull the image as follows:

```
VERSION="3.4.1"
docker pull genomicslab/nanocaller:${VERSION}
```

And then run it as follows:

```
VERSION="3.4.1"
docker run \
-v 'YOUR_INPUT_DIR':'/input/' \
-v 'YOUR_WORKING_DIR':'/output/' \
genomicslab/nanocaller:${VERSION} \
NanoCaller \
--bam /input/YOUR_BAM \
--ref /input/YOUR_REF \
--output /output \
--cpu NUMBER_OF_CPUS_TO_USE
--preset PRESET
```

## Running NanoCaller with Singularity

You can use docker image for NanoCaller with Singularity with a similar syntax as described for docker above. Please check the [NanoCaller Docker Hub repository](https://hub.docker.com/repository/docker/genomicslab/nanocaller) for the most up to date version of NanoCaller docker image. You can pull the image as follows:

```
VERSION="3.4.1"
singularity pull docker://genomicslab/nanocaller:${VERSION}
```

And then run the image as follows:
```
VERSION="3.4.1"
singularity run \
--bind 'YOUR_INPUT_DIR':'/input/' \
--bind 'YOUR_WORKING_DIR':'/output/' \
--pwd /app nanocaller_${VERSION}.sif \
--bam /input/YOUR_BAM \
--ref /input/YOUR_REF \
--output /output \
--cpu NUMBER_OF_CPUS_TO_USE
--preset PRESET
```

## Running NanoCaller with Bioconda Installation

NanoCaller can call variants from whole genome or several chromosomes using a single command `NanoCaller`. Assuming all your input files are in a folder `YOUR_INPUT_DIR`, and you want to use `YOUR_OUTPUT_DIR` to store the results. You can use `--regions` or `--bed` to specify regions for variant calling, otherwise all contigs in the input BAM file will be used.

```
NanoCaller \
--bam YOUR_BAM \
--ref YOUR_REF \
--output OUTPUT_DIRECTORY \
--cpu NUMBER_OF_CPUS_TO_USE
--preset PRESET
```

## General Usage Options

```
usage: NanoCaller [-h] [--mode {snps,indels,all}] [--sequencing {ont,ul_ont,ul_ont_extreme,pacbio}] [--cpu CPU] [--mincov MINCOV] [--maxcov MAXCOV]
                  [--suppress_progress_bar] [--output OUTPUT] [--prefix PREFIX] [--sample SAMPLE] [--regions [REGIONS [REGIONS ...]]] [--bed BED]
                  [--wgs_contigs {chr1-22XY,1-22XY}] [--exclude_bed {hg38,hg19,mm10,mm39}] [--preset {ont,ul_ont,ul_ont_extreme,ccs,clr}] --bam BAM --ref
                  REF [--snp_model SNP_MODEL] [--min_allele_freq MIN_ALLELE_FREQ] [--min_nbr_sites MIN_NBR_SITES] [--neighbor_threshold NEIGHBOR_THRESHOLD]
                  [--supplementary] [--indel_model INDEL_MODEL] [--ins_threshold INS_THRESHOLD] [--del_threshold DEL_THRESHOLD] [--win_size WIN_SIZE]
                  [--small_win_size SMALL_WIN_SIZE] [--impute_indel_phase] [--phase] [--enable_whatshap]

optional arguments:
  -h, --help            show this help message and exit

Required Arguments:
  --bam BAM             Bam file, should be phased if 'indel' mode is selected (default: None)
  --ref REF             Reference genome file with .fai index (default: None)

Preset:
  --preset {ont,ul_ont,ul_ont_extreme,ccs,clr}
                        Apply recommended preset values for SNP and Indel calling parameters. 'ont' works well for all types of ONT sequencing datasets.
                        However, use 'ul_ont' if you have several ultra-long ONT reads up to 100kbp long, and 'ul_ont_extreme' if you have several ultra-
                        long ONT reads up to 300kbp long. For PacBio CCS (HiFi) and CLR reads, use 'ccs'and 'clr' respectively. Presets are described in
                        detail here: github.com/WGLab/NanoCaller/blob/master/docs/Usage.md#preset-options. (default: None)

Configurations:
  --mode {snps,indels,all}
                        NanoCaller mode to run. 'snps' mode quits NanoCaller without using WhatsHap for phasing. In this mode, if you want NanoCaller to
                        phase SNPs and BAM files, use --phase argument additionally. (default: all)
  --sequencing {ont,ul_ont,ul_ont_extreme,pacbio}
                        Sequencing type. 'ont' works well for any type of ONT sequencing datasets. However, use 'ul_ont' if you have several ultra-long ONT
                        reads up to 100kbp long, and 'ul_ont_extreme' if you have several ultra-long ONT reads up to 300kbp long. For PacBio CCS (HiFi) and
                        CLR reads, use 'pacbio'. (default: ont)
  --cpu CPU             Number of CPUs to use. (default: 1)
  --mincov MINCOV       Minimum coverage to call a variant (default: 4)
  --maxcov MAXCOV       Maximum coverage of reads to use. If sequencing depth at a candidate site exceeds maxcov then reads are downsampled. (default: 160)
  --suppress_progress_bar
                        Do not show progress bar. (default: False)

Variant Calling Regions:
  Use only one of these options to specify regions for variant calling: --regions or --bed option or --wgs_contigs. If none is provided then whole genome
  variant calling is assumed and all contigs in the BAM file are used.

  --regions [REGIONS [REGIONS ...]]
                        A space/whitespace separated list of regions specified as "CONTIG_NAME" or "CONTIG_NAME:START-END". If you want to use
                        "CONTIG_NAME:START-END" format then specify both start and end coordinates. For example: chr3 chr6:28000000-35000000 chr22.
                        (default: None)
  --bed BED             A BED file specifying regions for variant calling. (default: None)
  --wgs_contigs {chr1-22XY,1-22XY}
                        Preset list of chromosomes to use for variant calling on human genomes. "chr1-22XY" option will assume human reference genome with
                        "chr" prefix present in the chromosome notation, and run NanoCaller on chr1 to chr22, chrX and chrY. "1-22XY" option will assume no
                        "chr" prefix is present in the chromosome notation and run NanoCaller on chromosomes 1-22, X and Y. (default: None)
  --exclude_bed {hg38,hg19,mm10,mm39}
                        Path to bgzipped and tabix indexed BED file containing intervals to ignore for variant calling. BED files of centromere and
                        telomere regions for the following genomes are included in NanoCaller: hg38, hg19, mm10 and mm39. (default: None)

SNP Calling:
  --snp_model SNP_MODEL
                        NanoCaller SNP model to be used (default: ONT-HG002)
  --min_allele_freq MIN_ALLELE_FREQ
                        minimum alternative allele frequency (default: 0.15)
  --min_nbr_sites MIN_NBR_SITES
                        minimum number of nbr sites (default: 1)
  --neighbor_threshold NEIGHBOR_THRESHOLD
                        SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended, for
                        PacBio CCS anc CLR reads '0.3,0.7' and '0.3,0.6' are recommended respectively (default: 0.4,0.6)
  --supplementary       Use supplementary reads, not fully supported at the moment (default: False)

Indel Calling:
  --indel_model INDEL_MODEL
                        NanoCaller indel model to be used (default: ONT-HG002)
  --ins_threshold INS_THRESHOLD
                        Insertion Threshold (default: 0.4)
  --del_threshold DEL_THRESHOLD
                        Deletion Threshold (default: 0.6)
  --win_size WIN_SIZE   Size of the sliding window in which the number of indels is counted to determine indel candidate site. Only indels longer than 2bp
                        are counted in this window. Larger window size can increase recall, but use a maximum of 50 only (default: 40)
  --small_win_size SMALL_WIN_SIZE
                        Size of the sliding window in which indel frequency is determined for small indels (default: 4)
  --impute_indel_phase  Infer read phase by rudimentary allele clustering if the no or insufficient phasing information is available, can be useful for
                        datasets without SNPs or regions with poor phasing quality (default: False)

Output Options:
  --output OUTPUT       VCF output path, default is current working directory (default: None)
  --prefix PREFIX       VCF file prefix (default: variant_calls)
  --sample SAMPLE       VCF file sample name (default: SAMPLE)

Phasing:
  --phase               Phase SNPs and BAM files if snps mode is selected. (default: False)
  --enable_whatshap     Allow WhatsHap to change SNP genotypes when phasing using --distrust-genotypes and --include-homozygous flags (this is not the same
                        as regenotyping), increasing the time needed for phasing. It has a negligible effect on SNP calling accuracy for Nanopore reads,
                        but may make a small improvement for PacBio reads. By default WhatsHap will only phase SNP calls produced by NanoCaller, but not
                        change their genotypes. (default: False)
```

## Understanding NanoCaller Output

Depending upon which mode is run, NanoCaller will produce the following files:

- PREFIX.snps.vcf.gz contains unphased SNP calls made by NanoCaller using a deep learning model. 
- PREFIX.snps.phased.vcf.gz contains SNP calls from PREFIX.snps.vcf.gz that are phased with WhatsHap if `all` mode is selected or `--phase` option was selected with `snps` mode. By default SNP calls in this file have the same genotype as in PREFIX.snps.vcf.gz file, unless `--enable_whatshap` flag is set which can allow WhatsHap to change genotypes.
- PREFIX.indels.vcf.gz contains indel calls made by NanoCaller using multiple sequence alignment. Some of these indel calls might be combined with nearby substitutions or multi-nucleotide substitutions.
- PREFIX.vcf.gz contains SNP calls from PREFIX.snps.phased.vcf.gz and indel calls from PREFIX.indels.vcf.gz.
- intermediate_phase_files directory will contain phased BAM files if `all` mode is selected or `--phase` option was selected with `snps` mode.

## Preset Options
Users can select recommended NanoCaller settings for various sequencing types by using preset option `--preset`. There are five presets, which are described below with their equivalent parameter settings. Any parameters included in the presets can be overwritten by specifying the parameter.

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
--del_threshold 0.4
--enable_whatshap
--impute_indel_phase
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
- `-nbr_t` option is sensitive to sequencing type so choose this accordingly.
- `-ins_t` and `del_t` are insertion and deletion frequency thresholds are per haplotype. Default values are slightly higher due to high error in ONT reads, but these thresholds can be lowered for CCS reads.
- CLR reads have incredibly low insertion and deletion freqencies in a pileup due highly variable placement of indels by aligners. To detect indels on CLR reads, you would need to set low frequency threhsolds, which leads to a huge increase in runtime.
- `--enable_whatshap` flag can improve SNP calling performance for PacBio reads.
