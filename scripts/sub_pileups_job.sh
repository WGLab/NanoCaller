#!/bin/bash



for i in {1..4}; do echo "bash /home/ahsanm1/repos/NanoVar/scripts/make_candidates.sh $i training 16" | qsub -q gpu.q -l GPU=0 -V -cwd  -pe smp 16  -N hg1tr.c$i -e /home/ahsanm1/dt_Nanovar/pileups/logs/hg1tr.c$i -o /home/ahsanm1/dt_Nanovar/pileups/logs/hg1tr.c$i;done


for i in {6..15}; do echo "bash /home/ahsanm1/repos/NanoVar/scripts/make_candidates.sh $i training 8" | qsub -q gpu.q -l GPU=0 -V -cwd -pe smp 8  -N hg1tr.c$i -e /home/ahsanm1/dt_Nanovar/pileups/logs/hg1tr.c$i -o /home/ahsanm1/dt_Nanovar/pileups/logs/hg1tr.c$i;done

