import pysam, sys, argparse

if __name__=='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--vcf')
    parser.add_argument('--bam')
    parser.add_argument('--output')
    
    args = parser.parse_args()
    
    vcf_path=args.vcf
    bam_path=args.bam
    output_path=args.output

    vcf = pysam.VariantFile(vcf_path)

    variant_alleles={}
    for rec in vcf.fetch():
        if max(len(x) for x in  rec.alleles )==1:
            variant_alleles[(rec.contig, rec.pos)]=rec.alleles

    bam=pysam.AlignmentFile(bam_path, 'rb')

    with open(output_path,'w') as outfile:
        for pcol in bam.pileup(min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800):
            rec_id=(pcol.reference_name, pcol.pos+1)
            if rec_id in variant_alleles:
                rec_alleles=variant_alleles[rec_id]
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()

                reads_per_var={x:[] for x in 'ACGT*'}

                for n,s in zip(name,seq):
                    reads_per_var[s].append(n)

                allele_str='\t'.join('{}:{}'.format(k,','.join(reads_per_var[k])) for k in rec_alleles)
                rec_str='{}\t{}\t{}\n'.format(rec_id[0], rec_id[1], allele_str)
                outfile.write(rec_str)