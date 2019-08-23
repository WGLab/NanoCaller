#!/bin/bash



bam=/home/ahsanm1/wglab_shared/datasets/NA12878_HG001_data/ont_rel6/alignments/GRCh38/HG001.rel6.hg38.sorted.bam
vcf=/home/ahsanm1/umair_wlab/data/ground_truth_SNPs/HG001.GRCh38.SNPonly.uniq.normalizeGT.vcf.gz
ref=/home/ahsanm1/umair_wlab/data/GRCh38.fa
bed=/home/ahsanm1/dt_Nanovar/bed_by_chrom/HG001/



if [ "$2" = "training" ];then

output=/home/ahsanm1/dt_Nanovar/pileups/HG001_ont/training/chr$1

python /home/ahsanm1/repos/NanoVar/scripts/generate_pileups.py -r chr$1 -m $2 -t 0.0 -bam $bam -ref $ref -vcf $vcf -o $output -w 16 -d 32 -cpu $3 -mincov 8 -bed $bed

else

output=/home/ahsanm1/dt_Nanovar/pileups/HG001_ont/testing/$4/chr$1

python /home/ahsanm1/repos/NanoVar/scripts/generate_pileups.py -r chr$1 -m $2 -t 0.$4 -bam $bam -ref $ref -vcf $vcf -o $output -w 16 -d 32 -cpu $3 -mincov 8

fi


echo '
HG1 ont
'
