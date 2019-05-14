#!/bin/bash

python CNN_model.py -r 0.001 -i 20 -s 20 -n 3 -w 20 -train /home/ahsanm1/projects/variant_calling/Pileup/output/chr20_48000000_53000000_pileups.gz -test /home/ahsanm1/projects/variant_calling/Pileup/output/chr20_43000000_48000000_pileups.gz -p 0 -model /home/ahsanm1/projects/variant_calling/Pileup/models/cnn_model -m train
