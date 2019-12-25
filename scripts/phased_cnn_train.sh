#!/bin/bash

genome=$1
model=/home/ahsanm/dt_Nanovar/models/haploid_models/${genome}_ont/

mkdir -p $model

python phased_model_run.py -r 0.001 -i 100 -s 20 -train /home/ahsanm/dt_Nanovar/pileups/${genome}_ont/phased/training -test /home/ahsanm/dt_Nanovar/pileups/HG2_ont/phased/training/chr22/chr22.pileups.:/home/ahsanm/dt_Nanovar/pileups/HG1_ont/phased/training/chr21/chr21.pileups. -model $model/model -m train -dim 32:9:7 -vcf None -cpu 8 -val 1 -rt 0  -chrom 2:22 -ratio 4 -neg  neg.combined -rt_path $model

