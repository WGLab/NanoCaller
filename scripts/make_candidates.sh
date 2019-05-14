#!/bin/bash

python generate_candidate_pileups.py -r chr20:$1-$2 -m training -bam /home/ahsanm1/projects/variant_calling/data/nanopore_chr20.bam -ref /home/ahsanm1/projects/variant_calling/data/ref_GRCh38_chr20.fa -vcf /home/ahsanm1/shared_data/NA12878/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -o /home/ahsanm1/projects/variant_calling/Pileup/output
