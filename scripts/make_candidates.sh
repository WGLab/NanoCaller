#!/bin/bash

python /home/ahsanm1/repos/NanoVar/scripts/generate_pileups.py -r chr$1 -m $2 -t 0.15 -bam /home/ahsanm1/shared_data/NA12878/Nanopore/Alignmentsbychromosome/chr$1.sorted.bam -ref /home/ahsanm1/projects/variant_calling/data/GRCh38.fa -vcf /home/ahsanm1/umair_wlab/data/NanoVar_data/HG001.inbed.SNPonly.uniq.normalizeGT.vcf.gz -o /home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/$2/chr$1 -w 16 -d 32 -cpu $3 -bed True


echo '
python /home/ahsanm1/repos/NanoVar/scripts/generate_pileups.py -r chr$1 -m $2 -t 0.15 -bam /home/ahsanm1/shared_data/NA12878/Nanopore/Alignmentsbychromosome/chr$1.sorted.bam -ref /home/ahsanm1/projects/variant_calling/data/GRCh38.fa -vcf /home/ahsanm1/umair_wlab/data/NanoVar_data/HG001.inbed.SNPonly.uniq.normalizeGT.vcf.gz -o /home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/$2/chr$1 -w 16 -d 32 -cpu $3 -bed True
'
