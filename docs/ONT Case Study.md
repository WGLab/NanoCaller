# ONT Case Study
## HG002
```
# Download reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz \
-O -| gunzip -c > GRCh38.fa

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai \
-O GRCh38.fa.fai

# run `pip install awscli` if you do not have aws cli API installed.
# download FASTQ files

for i in {1..3};do aws s3 cp s3://human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/GM24385_${i}_Guppy_4.2.2_prom.fastq.gz \
./ --no-sign-request;done

CPU=32  # number of threads or parallel jobs

# align FASTQ files to reference genome, write sorted alignments to a BAM file
minimap2 -a -z 600,200 -x map-ont GRCh38.fa GM24385_1_Guppy_4.2.2_prom.fastq.gz \
GM24385_2_Guppy_4.2.2_prom.fastq.gz  GM24385_3_Guppy_4.2.2_prom.fastq.gz -t $CPU \
-v1|samtools view -Shu |samtools sort -@ $CPU -o HG002.Guppy_4.2.2_prom.bam --output-fmt BAM

# create index for the BAM file
samtools index HG002.Guppy_4.2.2_prom.bam -@ $CPU

# run nanocaller
VERSION=0.3.0
docker run -it -v ${PWD}:'/mnt/'  umahsn/nanocaller:${VERSION} python NanoCaller_WGS.py \
-bam /mnt/HG002.Guppy_4.2.2_prom.bam -ref /mnt/GRCh38.fa -prefix HG002 \
-o /mnt/calls -cpu $CPU --exclude_bed hg38


# If you want to run NanoCaller without docker, run the following command `python PATH_TO_NANOCALLER_REPO/NanoCaller_WGS.py -bam HG002.Guppy_4.2.2_prom.bam -ref GRCh38.fa -prefix HG002 -o calls --exclude_bed hg38`


# run `conda install bedtools` to install bedtools to create BED files for variant calling evaluation in difficult-to-map genomic regions.


mkdir -p evaluation_files/difficult_regions/  evaluation_files/difficult_regions/homopolymers analysis

# Create SDF (sequence data file) for reference genome required for RTG vcfeval 
rtg RTG_MEM=4G format -f fasta GRCh38.fa -o evaluation_files/GRCh38.sdf

# Download GIAB benchmark variants for HG002
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -P evaluation_files/
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi -P evaluation_files/
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -P evaluation_files/

# Evaluate variant calls using RTG vcfeval in high-confidence regions
rtg  RTG_MEM=12G vcfeval -b evaluation_files/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
-c calls/HG002.final.vcf.gz -t evaluation_files/GRCh38.sdf -e evaluation_files/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
-Z -o analysis/whole_genome -f 'QUAL'

# Display SNP and indels performance in high-confidence regions. First two results show performance with thresholds that give same
# highest F1 scores, but with higher recall and precision respectively. Last result shows performance when no quality threshold is used 
for TYPE in {snp,non_snp};do echo "$TYPE performance"; grep 'f_measure' analysis/whole_genome/${TYPE}_roc.tsv; \
tail -n +8 analysis/whole_genome/${TYPE}_roc.tsv |awk -v max=0 '{if($8>max){want=$0; max=$8}}END{print want}';\
tail -n +8 analysis/whole_genome/${TYPE}_roc.tsv |awk -v max=0 '{if($8>=max){want=$0; max=$8}}END{print want}';\
tail -1 analysis/whole_genome/${TYPE}_roc.tsv; echo "";done


# Download BED files for difficult-to-map regions from GIAB genome stratification v2
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/union/GRCh38_alldifficultregions.bed.gz \
-O - |gzip -dc|grep -v ^'#' > evaluation_files/difficult_regions/GRCh38_alldifficultregions.bed

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/mappability/GRCh38_lowmappabilityall.bed.gz \
-O - |gzip -dc|grep -v ^'#' > evaluation_files/difficult_regions/GRCh38_lowmappabilityall.bed

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/SegmentalDuplications/GRCh38_segdups.bed.gz \
-O - |gzip -dc|grep -v ^'#' > evaluation_files/difficult_regions/GRCh38_segdups.bed

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz \
-O - |gzip -dc|grep -v ^'#' > evaluation_files/difficult_regions/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/OtherDifficult/GRCh38_MHC.bed.gz \
 -O - |gzip -dc|grep -v ^'#' > evaluation_files/difficult_regions/GRCh38_MHC.bed

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/LowComplexity/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed.gz \
-O - |gzip -dc|grep -v ^'#' > evaluation_files/difficult_regions/homopolymers/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/LowComplexity/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz \
-O - |gzip -dc|grep -v ^'#' > evaluation_files/difficult_regions/homopolymers/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed

# Intersect HG002 high-confidence regions with various difficult-to-map regions
for x in evaluation_files/difficult_regions/GRCh38_*bed;do d="${x##*/}";\
bedtools intersect -a evaluation_files/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
-b $x > evaluation_files/difficult_regions/HG002_$d; done

# Remove homopolymer regions from HG002 high-confidence regions
bedtools subtract -a evaluation_files/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
-b evaluation_files/difficult_regions/homopolymers/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed| bedtools subtract -a - \
evaluation_files/difficult_regions/homopolymers/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed > \
evaluation_files/difficult_regions/HG002_minus_homopolymer_repeats.bed


# Evaluate performance in various difficult-to-map regions, and outside homopolymer repeats.
for x in evaluation_files/difficult_regions/HG002_*bed;do f="${x##*/}";d="${f%.bed}"; \
rtg  RTG_MEM=12G vcfeval -b evaluation_files/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
-c calls/HG002.final.vcf.gz -t evaluation_files/GRCh38.sdf -e $x -Z -o analysis/$d -f 'QUAL';done

# Show SNP and indels performance for whole genome, various difficult-to-map regions, and outside homopolymer repeats
# First two results show performance with highest F1 scores. Third result shows performance with quality score cut-off
# specified in variable THRESHOLD
TYPE=snp
THRESHOLD=162
for x in analysis/*;do f="${x##*/}"; echo "$f: $TYPE performance"; grep 'f_measure' $x/${TYPE}_roc.tsv; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>=max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=$THRESHOLD '{if($1>=max){want=$0}}END{print want}'; \
tail -1 $x/${TYPE}_roc.tsv; echo "";done

TYPE=non_snp
THRESHOLD=144
for x in analysis/*;do f="${x##*/}"; echo "$f: $TYPE performance"; grep 'f_measure' $x/${TYPE}_roc.tsv; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>=max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=$THRESHOLD '{if($1>=max){want=$0}}END{print want}'; \
tail -1 $x/${TYPE}_roc.tsv; echo "";done
```
