# ONT Variant Calling

## Alignment Commands

mnimap2 -a -z 600,200 -x map-ont GRCh38.fa HG002_GM24385_1_2_3_Guppy_3.6.0_prom.fastq.gz

## NanoCaller Commands

NanoCaller v0.1.0 (commit aa62f3b8336468767808f51807b129975f6c7924): https://github.com/WGLab/NanoCaller/releases/tag/v0.1.0

`for i in {1..22};do python NanoCaller.py -chrom chr$i -bam HG002_ont.bam -ref GRCh38.fa -prefix HG002.chr$i -mode both -seq ont -nbr_t '0.4,0.6' -ins_t 0.4 -del_t 0.6 -model NanoCaller2 -vcf chr$i -cpu 16;done`

`ls -1 chr*/*phased.vcf.gz|bcftools concat -f - -a| bcftools sort| bgziptabix HG002_nanocaller.snps.vcf.gz`

`ls -1 chr*/*.decomposed.vcf.gz|bcftools concat -f - -a|rtg vcffilter -i - --non-snps-only -o -| bcftools sort| bgziptabix HG002_nanocaller.indels.vcf.gz`


Min-max normalize scores and remove outlier scores

`python nanocaller_snp_vcf_fixer.py HG002_nanocaller.snps.vcf.gz nanocaller |bgziptabix  HG002_nanocaller.fixed.snps.vcf.gz`

`python nanocaller_indel_fixer.py HG002_nanocaller.indels.vcf.gz nanocaller |bgziptabix  HG002_nanocaller.fixed.indels.vcf.gz`


## Clair Commands

Clair v2.0.1 (commit f97b099d00188dfdcdeb4bcc32080e1a8f457d5d): https://github.com/HKU-BAL/Clair/releases/tag/v2.0.1
Clair model used: http://www.bio8.cs.hku.hk/clair_models/ont/1_124x.tar

`pypy3 clair.py callVarBamParallel \
       --chkpnt_fn model \
       --ref_fn GRCh38.fa \
       --bam_fn HG002_ont.bam \
       --sampleName HG002 \
       --output_prefix HG002 \
       > commands.sh`
       
cat command.sh | parallel -j16

vcfcat HG002*.vcf|bcftools sort|bgziptabix HG002_clair.vcf.gz

rtg vcffilter -i HG002_clair.vcf.gz --snps-only -o -|bgziptabix HG002_clair.snps.vcf.gz
rtg vcffilter -i HG002_clair.vcf.gz --non-snps-only -o -|bgziptabix HG002_clair.indels.vcf.gz

python clair_snp_vcf_fixer.py  HG002_clair.snps.vcf.gz clair |bgziptabix  HG002_clair.fixed.snps.vcf.gz
python clair_indels_vcf_fixer.py  HG002_clair.indels.vcf.gz clair |bgziptabix  HG002_clair.fixed.indels.vcf.gz

## Medaka Commands

Medaka v0.10.0 (commit eeb1af187c1aafc693ea19b25795485d385461eb): https://github.com/nanoporetech/medaka/releases/tag/v0.10.0


for i in {1..22};do mkdir chr$i; medaka_variant -i HG002_ont.bam -r chr$i -f /home/ahsanm1/umair_wlab/data/GRCh38.fa -o chr$i -t 16 -d;done

ls -1 chr*/round_1_phased.vcf|bcftools concat -f - -a|bcftools sort| rtg vcfdecompose -t GRCh38.sdf -i - -o - --break-mnps --break-indels| bgziptabix HG002_medaka.temp.vcf.gz

python remove_bad_medaka_entries.py HG002_medaka.temp.vcf.gz |bgziptabix HG002_medaka.temp.vcf.gz

rtg vcffilter -i HG002_medaka.vcf.gz --snps-only -o -|bgziptabix HG002_medaka.snps.vcf.gz
rtg vcffilter -i HG002_medaka.vcf.gz --non-snps-only -o -|bgziptabix HG002_medaka.indels.vcf.gz

python medaka_snp_vcf_fixer.py  HG002_medaka.snps.vcf.gz medaka |bgziptabix  HG002_medaka.fixed.snps.vcf.gz
python medaka_indels_vcf_fixer.py  HG002_medaka.indels.vcf.gz medaka |bgziptabix  HG002_medaka.fixed.indels.vcf.gz


## Ensemble Calling


bcftools merge HG002_nanocaller.fixed.snps.vcf.gz HG002_clair.fixed.snps.vcf.gz HG002_medaka.fixed.snps.vcf.gz | bgziptabix HG002_merged.snps.vcf.gz

python snp_ensemble.py HG002_merged.snps.vcf.gz ensemble | bzgiptabix HG002_ensemble.snps.vcf.gz

bcftools view HG002_ensemble.snps.vcf.gz -i "QUAL>=7.41" | bzgiptabix HG002_ensemble.snps.filtered.vcf.gz 

For HG003 and HG004 we used 15.6 and 13.6 as QUAL cut off for SNPs in the step above.


bcftools merge HG002_nanocaller.fixed.indels.vcf.gz HG002_clair.fixed.indels.vcf.gz HG002_medaka.fixed.indels.vcf.gz | bgziptabix HG002_merged.indels.vcf.gz

python indel_ensemble.py HG002_merged.indels.vcf.gz ensemble | bzgiptabix HG002_ensemble.indels.vcf.gz

### Final VCF file
bcftools concat -a HG002_ensemble.snps.filtered.vcf.gz HG002_ensemble.indels.vcf.gz| bgziptabix HG2_ont.vcf.gz


# PacBio Variant Calling


## Alignment Commands
Minimap2 parameters: -a -x map-pb -k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k --eqx --secondary=no GRCh38.fa HG002_35x_PacBio_14kb-15kb.fastq.gz



## NanoCaller Commands


`for i in {1..22};do python NanoCaller.py -chrom chr$i -bam HG002_CCS.bam -ref GRCh38.fa -prefix HG002.chr$i -mode both -seq pacbio -nbr_t '0.3,0.7' -ins_t 0.3 -del_t 0.3 -model NanoCaller1 -vcf chr$i -cpu 16;done`

### Final VCF file
`ls -1 chr*/*final.vcf.gz|bcftools concat -f - -a| bcftools sort| bgziptabix HG2_CCS.vcf.gz`








