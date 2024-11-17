# ONT Case Study
## HG002

```
# create the working directory

mkdir NanoCaller_ONT_Case_Study
cd NanoCaller_ONT_Case_Study

# Clone NanoCaller repository, install all the packages needed for the case study

git clone https://github.com/WGLab/NanoCaller.git
conda env create -f NanoCaller/environment.yml
conda activate nanocaller_env
pip install awscli
conda install -y -c bioconda minimap2 bedtools 

# Download reference genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz \
-O -| gunzip -c > GRCh38.fa

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai \
-O GRCh38.fa.fai

# run `pip install awscli` if you do not have aws cli API installed.
# download FASTQ files

for i in {1..3};do aws s3 cp s3://human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/GM24385_${i}_Guppy_4.2.2_prom.fastq.gz \
./ --no-sign-request;done

CPU=16  # number of threads or parallel jobs. Please set this number higher or lower according to the number of processes allowed in your system.

# align FASTQ files to reference genome, write sorted alignments to a BAM file
minimap2 -a -z 600,200 -x map-ont GRCh38.fa GM24385_1_Guppy_4.2.2_prom.fastq.gz \
GM24385_2_Guppy_4.2.2_prom.fastq.gz  GM24385_3_Guppy_4.2.2_prom.fastq.gz -t $CPU \
-v1|samtools view -Shu |samtools sort -@ $CPU -o HG002.Guppy_4.2.2_prom.bam --output-fmt BAM

# create index for the BAM file
samtools index HG002.Guppy_4.2.2_prom.bam -@ $CPU

# run nanocaller
VERSION=3.4.1
docker run -it -v ${PWD}:'/mnt/'  genomicslab/nanocaller:${VERSION} NanoCaller \
--bam /mnt/HG002.Guppy_4.2.2_prom.bam --ref /mnt/GRCh38.fa --prefix HG002 --preset ont \
--output /mnt/calls --cpu $CPU --exclude_bed hg38 --wgs_contigs chr1-22XY


# If you want to run NanoCaller without docker, run the following command `NanoCaller --bam HG002.Guppy_4.2.2_prom.bam --ref GRCh38.fa --prefix HG002 --preset ont --output calls --exclude_bed hg38 --cpu $CPU --wgs_contigs chr1-22XY`


# run `conda install -c bioconda bedtools` to install bedtools to create BED files for variant calling evaluation in difficult-to-map genomic regions.


mkdir -p evaluation_files/difficult_regions/  evaluation_files/difficult_regions/homopolymers analysis

# Create SDF (sequence data file) for reference genome required for RTG vcfeval 
rtg RTG_MEM=4G format -f fasta GRCh38.fa -o evaluation_files/GRCh38.sdf
header="#score true_positives_baseline false_positives true_positives_call false_negatives precision sensitivity f_measure"

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
for TYPE in {snp,non_snp};do echo "$TYPE performance"; echo "$header"; \
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
-b evaluation_files/difficult_regions/homopolymers/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed > \
evaluation_files/difficult_regions/HG002_minus_homopolymer_repeats.bed


# Evaluate performance in various difficult-to-map regions, and outside homopolymer repeats.
for x in evaluation_files/difficult_regions/HG002_*bed;do f="${x##*/}";d="${f%.bed}"; \
rtg  RTG_MEM=12G vcfeval -b evaluation_files/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
-c calls/HG002.final.vcf.gz -t evaluation_files/GRCh38.sdf -e $x -Z -o analysis/$d -f 'QUAL';done

# Show SNP and indels performance for whole genome, various difficult-to-map regions, and outside homopolymer repeats
# First two results show performance with highest F1 scores. Third result shows performance with quality score cut-off
# specified in variable THRESHOLD
TYPE=snp
THRESHOLD=50
for x in analysis/*;do f="${x##*/}"; echo "$f: $TYPE performance"; echo "$header"; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>=max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=$THRESHOLD '{if($1>=max){want=$0}}END{print want}'; \
tail -1 $x/${TYPE}_roc.tsv; echo "";done

TYPE=non_snp
THRESHOLD=125
for x in analysis/*;do f="${x##*/}"; echo "$f: $TYPE performance"; echo "$header"; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>=max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=$THRESHOLD '{if($1>=max){want=$0}}END{print want}'; \
tail -1 $x/${TYPE}_roc.tsv; echo "";done
```

At the end, you should be able to get the following performances:

### SNPs
##### HG002_GRCh38_alldifficultregions: snp performance
| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
56.790|598306.45|17955.00|598602.00|29644.55|0.9709|0.9528|0.9618|
45.708|599730.98|19493.00|600025.00|28220.02|0.9685|0.9551|0.9618|
50.001|599198.03|18834.00|599492.00|28752.97|0.9695|0.9542|0.9618|
30.104|601709.00|22825.00|602002.00|26242.00|0.9635|0.9582|0.9608|

##### HG002_GRCh38_AllTandemRepeatsandHomopolymers_slop5: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
47.017|161475.60|8801.00|161727.00|21425.40|0.9484|0.8829|0.9145|
39.404|162086.77|9528.00|162337.00|20814.23|0.9446|0.8862|0.9145|
50.001|161242.74|8560.00|161494.00|21658.26|0.9497|0.8816|0.9144|
30.104|162797.00|10812.00|163048.00|20104.00|0.9378|0.8901|0.9133|

##### HG002_GRCh38_lowmappabilityall: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
68.717|188162.78|6028.00|188260.00|4477.22|0.9690|0.9768|0.9729|
63.907|188302.79|6176.00|188400.00|4337.21|0.9683|0.9775|0.9729|
50.009|188645.47|6666.00|188743.00|3994.53|0.9659|0.9793|0.9725|
30.122|189169.00|7762.00|189267.00|3471.00|0.9606|0.9820|0.9712|

##### HG002_GRCh38_MHC: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
44.278|19870.98|162.00|19845.00|294.02|0.9919|0.9854|0.9887|
44.278|19870.98|162.00|19845.00|294.02|0.9919|0.9854|0.9887|
50.358|19861.97|156.00|19836.00|303.03|0.9922|0.9850|0.9886|
30.151|19883.00|190.00|19857.00|282.00|0.9905|0.9860|0.9883|

##### HG002_GRCh38_segdups: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
100.429|116132.91|7345.00|116147.00|4815.09|0.9405|0.9602|0.9503|
84.339|116620.93|7885.00|116635.00|4327.07|0.9367|0.9642|0.9503|
50.009|117522.98|9209.00|117537.00|3425.02|0.9273|0.9717|0.9490|
30.111|117956.00|10426.00|117970.00|2992.00|0.9188|0.9753|0.9462|

##### HG002_minus_homopolymer_repeats: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
65.891|2438152.23|16419.00|2439245.00|24647.77|0.9933|0.9900|0.9917|
39.249|2443188.19|21542.00|2444282.00|19611.81|0.9913|0.9920|0.9917|
50.005|2441180.35|19059.00|2442275.00|21619.65|0.9923|0.9912|0.9917|
30.104|2444942.00|24731.00|2446035.00|17858.00|0.9900|0.9927|0.9914|

##### whole_genome: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
52.149|3306392.90|31433.00|3306540.00|58722.10|0.9906|0.9825|0.9866|
43.392|3309458.40|34583.00|3309604.00|55656.60|0.9897|0.9835|0.9866|
50.001|3307150.89|32117.00|3307298.00|57964.11|0.9904|0.9828|0.9866|
30.104|3314038.00|42069.00|3314184.00|51077.00|0.9875|0.9848|0.9861|

### Indels

##### HG002_GRCh38_alldifficultregions: non_snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
132.500|103066.52|67226.00|102820.00|265594.48|0.6047|0.2796|0.3824|
124.010|104146.69|71794.00|103903.00|264514.31|0.5914|0.2825|0.3824|
125.000|104048.65|71326.00|103805.00|264612.35|0.5927|0.2822|0.3824|
2.230|123363.00|1516637.00|123080.00|245298.00|0.0751|0.3346|0.1226|

##### HG002_GRCh38_AllTandemRepeatsandHomopolymers_slop5: non_snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
125.210|80901.33|64722.00|80717.00|257979.67|0.5550|0.2387|0.3339|
124.200|80988.42|65158.00|80804.00|257892.58|0.5536|0.2390|0.3339|
125.000|80926.35|64841.00|80742.00|257954.65|0.5546|0.2388|0.3339|
2.230|98176.00|1323152.00|97970.00|240705.00|0.0689|0.2897|0.1114|

##### HG002_GRCh38_lowmappabilityall: non_snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
132.660|6278.22|2376.00|6252.00|4128.78|0.7246|0.6033|0.6584|
132.650|6278.22|2377.00|6252.00|4128.78|0.7245|0.6033|0.6584|
125.050|6340.13|2541.00|6314.00|4066.87|0.7130|0.6092|0.6571|
2.230|7297.00|76109.00|7259.00|3110.00|0.0871|0.7012|0.1549|

##### HG002_GRCh38_MHC: non_snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
113.110|1054.95|480.00|1056.00|631.05|0.6875|0.6257|0.6552|
113.110|1054.95|480.00|1056.00|631.05|0.6875|0.6257|0.6552|
125.520|1034.48|453.00|1037.00|651.53|0.6960|0.6136|0.6522|
2.290|1131.00|4747.00|1132.00|555.00|0.1925|0.6708|0.2992|

##### HG002_GRCh38_segdups: non_snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
137.600|6476.15|2827.00|6435.00|4340.85|0.6948|0.5987|0.6432|
134.090|6507.12|2892.00|6465.00|4309.88|0.6909|0.6016|0.6432|
125.020|6576.06|3079.00|6534.00|4240.94|0.6797|0.6079|0.6418|
2.230|7244.00|76223.00|7189.00|3573.00|0.0862|0.6697|0.1527|

##### HG002_minus_homopolymer_repeats: non_snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
129.020|165029.36|31486.00|164828.00|42278.64|0.8396|0.7961|0.8173|
124.870|165477.97|32137.00|165276.00|41830.03|0.8372|0.7982|0.8173|
125.000|165462.96|32112.00|165261.00|41845.04|0.8373|0.7982|0.8173|
2.230|170042.00|237129.00|169816.00|37266.00|0.4173|0.8202|0.5532|

##### whole_genome: non_snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
133.320|229350.20|97642.00|229064.00|296058.80|0.7011|0.4365|0.5381|
131.940|229742.27|98704.00|229454.00|295666.73|0.6992|0.4373|0.5381|
125.000|231653.78|104668.00|231368.00|293755.22|0.6885|0.4409|0.5376|
2.230|261097.00|2733489.00|260756.00|264312.00|0.0871|0.4969|0.1482|
