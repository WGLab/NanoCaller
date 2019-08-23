#!/bin/bash


module load cudnn cuda91 BCFtools

export CUDA_VISIBLE_DEVICES=-1

chrom=chr16
j=20
name=hg2.$chrom.$j.22aug.32dp.vcf

model_path=/home/ahsanm1/dt_Nanovar/models/22Aug/2/model

python model_run.py -r 0.0 -i 0 -s 0 -train /home -test /home/ahsanm1/dt_Nanovar/pileups/old/HG002_ont/testing/$j/$chrom/$chrom.pileups.test -model $model_path -m test -dim 32:33:5 -vcf /home/ahsanm1/dt_Nanovar/vcf_files/$name -chrom $chrom -cpu 8 

bgzip /home/ahsanm1/dt_Nanovar/vcf_files/$name -f
bcftools index /home/ahsanm1/dt_Nanovar/vcf_files/$name.gz
