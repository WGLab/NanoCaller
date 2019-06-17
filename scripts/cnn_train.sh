#!/bin/bash

module load cudnn cuda91

python CNN_model.py -r 0.001 -i 400 -s 100 -train /home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/training/chr22/chr22 -test /home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/training/chr22/chr22- -model /home/ahsanm1/umair_wlab/data/NanoVar_data/models/june16/model -m train -dim 32:33:5 -vcf /home/ahsanm1/umair_wlab/data/NanoVar_data/ -cpu 30

echo '
python CNN_model.py -r 0.001 -i 30 -s 100 -train - -test /home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/training/chr22/chr22- -model /home/ahsanm1/umair_wlab/data/NanoVar_data/models/june16/model -m train -dim 32:33:5 -vcf /home/ahsanm1/umair_wlab/data/NanoVar_data/ -cpu 28

'

echo 'checking new pileups'
