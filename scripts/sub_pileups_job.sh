#!/bin/bash

for i in {1..18}; do echo "bash /home/ahsanm1/repos/NanoVar/scripts/make_candidates.sh $i training 10"| qsub -V -cwd -l h_vmem=16G -pe smp 10  -N tr_chr$i -e tr_chr$i.o -o tr_chr$i.o;done

for i in {1..18}; do echo "bash /home/ahsanm1/repos/NanoVar/scripts/make_candidates.sh $i testing 10"| qsub -V -cwd -l h_vmem=16G -pe smp 10  -N test_chr$i -e test_chr$i.o -o test_chr$i.o;done


