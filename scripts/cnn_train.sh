#!/bin/bash

module load cudnn cuda91

export CUDA_VISIBLE_DEVICES=1
 
python model_run-Copy1.py -chrom 18 -r 0.001 -i 60 -s 50 -train /home/ahsanm1/umair_wlab/data/NA24385_nanopore/pacbio/pileups/training -test /home/ahsanm1/umair_wlab/data/NA24385_nanopore/pileups/training/chr22/chr22_pileups_ -model /home/ahsanm1/umair_wlab/data/NA24385_nanopore/models/july15/pcb/model -m train -dim 32:33:5 -vcf None -cpu 12 -val 0
