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
THRESHOLD=162
for x in analysis/*;do f="${x##*/}"; echo "$f: $TYPE performance"; echo "$header"; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=0 '{if($8>=max){want=$0; max=$8}}END{print want}'; \
tail -n +8 $x/${TYPE}_roc.tsv |awk -v max=$THRESHOLD '{if($1>=max){want=$0}}END{print want}'; \
tail -1 $x/${TYPE}_roc.tsv; echo "";done

TYPE=non_snp
THRESHOLD=144
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
| 156.333        | 596478.27   | 19602  | 596715  | 31472.73 | 0.9682    | 0.9499 | 0.959  |
| 142.715        | 597792.69   | 21029  | 598029  | 30158.31 | 0.966     | 0.952  | 0.959  |
| 162.003        | 595897.28   | 19074  | 596135  | 32053.72 | 0.969     | 0.949  | 0.9589 |
| 30.105         | 604957      | 43263  | 605186  | 22994    | 0.9333    | 0.9634 | 0.9481 |


##### HG002_GRCh38_AllTandemRepeatsandHomopolymers_slop5: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 120.972        | 161087.05   | 9925   | 161283  | 21813.95 | 0.942     | 0.8807 | 0.9104 |
| 117.115        | 161258.98   | 10131  | 161455  | 21642.02 | 0.941     | 0.8817 | 0.9104 |
| 162.003        | 158846.5    | 8017   | 159042  | 24054.5  | 0.952     | 0.8685 | 0.9083 |
| 30.106         | 164605      | 18423  | 164792  | 18296    | 0.8994    | 0.9    | 0.8997 |



##### HG002_GRCh38_lowmappabilityall: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 181.082        | 187662.88   | 5976   | 187760  | 4977.12  | 0.9692    | 0.9742 | 0.9717 |
| 165.09         | 188088.4    | 6427   | 188186  | 4551.6   | 0.967     | 0.9764 | 0.9717 |
| 162.044        | 188162.07   | 6519   | 188260  | 4477.93  | 0.9665    | 0.9768 | 0.9716 |
| 30.105         | 189943      | 13381  | 190042  | 2697     | 0.9342    | 0.986  | 0.9594 |


##### HG002_GRCh38_MHC: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 154.099        | 19818.84    | 189    | 19783   | 346.16   | 0.9905    | 0.9828 | 0.9867 |
| 154.099        | 19818.84    | 189    | 19783   | 346.16   | 0.9905    | 0.9828 | 0.9867 |
| 162.596        | 19805.81    | 183    | 19770   | 359.19   | 0.9908    | 0.9822 | 0.9865 |
| 30.393         | 19904       | 451    | 19868   | 261      | 0.9778    | 0.9871 | 0.9824 |


##### HG002_GRCh38_segdups: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 184.36         | 116436.73   | 7888   | 116451  | 4511.27  | 0.9366    | 0.9627 | 0.9495 |
| 183.506        | 116467.73   | 7922   | 116482  | 4480.27  | 0.9363    | 0.963  | 0.9495 |
| 162.02         | 116968.76   | 8626   | 116983  | 3979.24  | 0.9313    | 0.9671 | 0.9489 |
| 30.105         | 118581      | 15725  | 118595  | 2367     | 0.8829    | 0.9804 | 0.9291 |


##### HG002_minus_homopolymer_repeats: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 208.453        | 2430326.94  | 27820  | 2431360 | 32473.06 | 0.9887    | 0.9868 | 0.9878 |
| 192.717        | 2433899.16  | 31482  | 2434935 | 28900.84 | 0.9872    | 0.9883 | 0.9878 |
| 162            | 2439048.52  | 39384  | 2440084 | 23751.48 | 0.9841    | 0.9904 | 0.9872 |
| 30.106         | 2450791     | 114883 | 2451831 | 12009    | 0.9552    | 0.9951 | 0.9748 |


##### whole_genome: snp performance

| Quality_Cutoff | TP_baseline | FP     | TP_call | FN       | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 183.284        | 3294583.7   | 45810  | 3294674 | 70531.3  | 0.9863    | 0.979  | 0.9827 |
| 171.709        | 3298466.2   | 49830  | 3298554 | 66648.8  | 0.9851    | 0.9802 | 0.9827 |
| 162            | 3301265.86  | 53468  | 3301353 | 63849.14 | 0.9841    | 0.981  | 0.9825 |
| 30.104         | 3324244     | 156174 | 3324321 | 40871    | 0.9551    | 0.9879 | 0.9712 |


### Indels
##### HG002_GRCh38_alldifficultregions: non_snp performance

| Quality_Cutoff | TP_baseline | FP      | TP_call | FN        | Precision | Recall | F1     |
|----------------|-------------|---------|---------|-----------|-----------|--------|--------|
| 260.57         | 86646.32    | 125755  | 86274   | 282014.68 | 0.4069    | 0.235  | 0.298  |
| 259.53         | 86997.38    | 127760  | 86624   | 281663.62 | 0.4041    | 0.236  | 0.298  |
| 144            | 108498.34   | 371321  | 108048  | 260162.66 | 0.2254    | 0.2943 | 0.2553 |
| 2.23           | 115645      | 860882  | 115154  | 253016    | 0.118     | 0.3137 | 0.1715 |

##### HG002_GRCh38_AllTandemRepeatsandHomopolymers_slop5: non_snp performance

| Quality_Cutoff | TP_baseline | FP      | TP_call | FN        | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 257.23         | 65564.57    | 127868  | 65254   | 273316.43 | 0.3379    | 0.1935 | 0.2461 |
| 249.16         | 67828.75    | 143930  | 67506   | 271052.25 | 0.3193    | 0.2002 | 0.2461 |
| 144            | 83918.67    | 359887  | 83545   | 254962.33 | 0.1884    | 0.2476 | 0.214  |
| 2.23           | 90283       | 821504  | 89885   | 248598    | 0.0986    | 0.2664 | 0.144  |

##### HG002_GRCh38_lowmappabilityall: non_snp performance

| Quality_Cutoff | TP_baseline | FP      | TP_call | FN        | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 252.8          | 6032.4      | 2876    | 6003    | 4374.6    | 0.6761    | 0.5796 | 0.6242 |
| 248.44         | 6100.31     | 3024    | 6069    | 4306.69   | 0.6674    | 0.5862 | 0.6242 |
| 144.04         | 6866.37     | 7673    | 6828    | 3540.63   | 0.4709    | 0.6598 | 0.5495 |
| 2.23           | 7242        | 24006   | 7195    | 3165      | 0.2306    | 0.6959 | 0.3464 |

##### HG002_GRCh38_MHC: non_snp performance

| Quality_Cutoff | TP_baseline | FP      | TP_call | FN        | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 215.99         | 990.37      | 363     | 999     | 695.63    | 0.7335    | 0.5874 | 0.6524 |
| 215.99         | 990.37      | 363     | 999     | 695.63    | 0.7335    | 0.5874 | 0.6524 |
| 144.25         | 1051.24     | 575     | 1060    | 634.76    | 0.6483    | 0.6235 | 0.6357 |
| 2.23           | 1101        | 1530    | 1106    | 585       | 0.4196    | 0.653  | 0.5109 |

##### HG002_GRCh38_segdups: non_snp performance

| Quality_Cutoff | TP_baseline | FP      | TP_call | FN        | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 258.97         | 5936.1      | 3115    | 5886    | 4880.9    | 0.6539    | 0.5488 | 0.5968 |
| 258.89         | 5937.1      | 3117    | 5887    | 4879.9    | 0.6538    | 0.5489 | 0.5968 |
| 144.03         | 6812.58     | 8266    | 6752    | 4004.42   | 0.4496    | 0.6298 | 0.5247 |
| 2.23           | 7121        | 24005   | 7051    | 3696      | 0.227     | 0.6583 | 0.3376 |

##### HG002_minus_homopolymer_repeats: non_snp performance

| Quality_Cutoff | TP_baseline | FP      | TP_call | FN        | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 143.45         | 159931.17   | 33199   | 159624  | 47376.83  | 0.8278    | 0.7715 | 0.7987 |
| 125.04         | 160660.28   | 34296   | 160355  | 46647.72  | 0.8238    | 0.775  | 0.7987 |
| 144.01         | 159903.17   | 33167   | 159596  | 47404.83  | 0.8279    | 0.7713 | 0.7986 |
| 2.23           | 163620      | 48092   | 163301  | 43688     | 0.7725    | 0.7893 | 0.7808 |

##### whole_genome: non_snp performance

| Quality_Cutoff | TP_baseline | FP      | TP_call | FN        | Precision | Recall | F1     |
|----------------|-------------|--------|---------|----------|-----------|--------|--------|
| 260.49         | 207019.9    | 146415  | 206586  | 318389.1  | 0.5852    | 0.394  | 0.471  |
| 259.86         | 207441.91   | 147783  | 207006  | 317967.09 | 0.5835    | 0.3948 | 0.471  |
| 144            | 243167.44   | 432410  | 242643  | 282241.56 | 0.3594    | 0.4628 | 0.4046 |
| 2.23           | 254059      | 1081248 | 253489  | 271350    | 0.1899    | 0.4835 | 0.2727 |
