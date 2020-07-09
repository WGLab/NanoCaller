import time, argparse, os

if __name__ == '__main__':
    t=time.time()
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('Required arguments')
    parser.add_argument("-mode",  "--mode",  help="Testing mode, options are 'snps', 'indels' and 'both'", type=str, default='both')
    parser.add_argument("-seq",  "--sequencing",  help="Sequencing type, options are 'ont' and 'pacbio'", type=str, default='ont')
    parser.add_argument("-model",  "--model",  help="NanoCaller SNP model to be used, options are 'NanoCaller1' (trained on HG001 Nanopore reads), 'NanoCaller2' (trained on HG002 Nanopore reads) and 'NanoCaller3' (trained on HG003 PacBio reads) ", default='NanoCaller1')
    parser.add_argument("-vcf",  "--vcf",  help="VCF output path", type=str, default='')
    
    requiredNamed.add_argument("-chrom",  "--chrom",  help="Chromosome", required=True)
    parser.add_argument("-cpu",  "--cpu",  help="CPUs", type=int, default=1)
    parser.add_argument("-min_allele_freq",  "--min_allele_freq",  help="minimum alternative allele frequency", type=float,  default=0.15)
    parser.add_argument("-min_nbr_sites",  "--min_nbr_sites",  help="minimum number of nbr sites", type=int,  default =1)
    
    requiredNamed.add_argument("-bam",  "--bam",  help="Bam file", required=True)
    requiredNamed.add_argument("-ref",  "--ref",  help="reference genome file with .fai index", required=True)
    
    requiredNamed.add_argument("-prefix",  "--prefix",  help="VCF file prefix", type=str, required=True)
    parser.add_argument("-sample",  "--sample",  help="VCF file sample name", type=str, default='SAMPLE')
    
    parser.add_argument("-mincov",  "--mincov",  help="min coverage", type=int, default=8)
    parser.add_argument("-maxcov",  "--maxcov",  help="max coverage", type=int, default=160)
    
    parser.add_argument("-start",  "--start",  help="start, default is 1", type=int)
    parser.add_argument("-end",  "--end",  help="end, default is the end of contig", type=int)
    parser.add_argument("-nbr_t",  "--neighbor_threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended and for PacBio reads '0.3,0.7' is recommended", type=str, default='0.4,0.6')
    parser.add_argument("-ins_t", "--ins_threshold", help="Insertion Threshold",type=float,default=0.4)
    parser.add_argument("-del_t", "--del_threshold", help="Deletion Threshold",type=float,default=0.6)
    
    args = parser.parse_args()
    
    
    
    try:  
        os.mkdir(args.vcf)  
    except OSError as error:  
        pass 
    
    
    with open(os.path.join(args.vcf,'args'),'w') as file:
        file.write(str(args))
        
    import snpCaller, indelCaller
    if not args.end:
        try:
            with open(args.ref+'.fai','r') as file:
                for line in file:
                    if line.split('\t')[0]==args.chrom:
                        end=int(line.split('\t')[1])

        except FileNotFoundError:
            print('index file .fai required for reference genome file', flush=True)
            sys.exit()
            
    else:
        end=args.end
    
    if not args.start:
        start=1
    else:
        start=args.start
    
    threshold=[float(args.neighbor_threshold.split(',')[0]), float(args.neighbor_threshold.split(',')[1])]
    
    in_dict={'chrom':args.chrom, 'start':start, 'end':end, 'sam_path':args.bam, 'fasta_path':args.ref, \
             'mincov':args.mincov,  'maxcov':args.maxcov, 'min_allele_freq':args.min_allele_freq, 'min_nbr_sites':args.min_nbr_sites, \
             'threshold':threshold, 'model':args.model, 'cpu':args.cpu,  'vcf_path':args.vcf,'prefix':args.prefix,'sample':args.sample, \
            'seq':args.sequencing}
        
    snp_vcf=''
    if args.mode in ['snps','both']:
        
        snp_vcf=snpCaller.test_model(in_dict)

        if snp_vcf:
            stream=os.popen("whatshap phase %s.vcf.gz %s -o %s.phased.preclean.vcf -r %s --ignore-read-groups --chromosome %s --distrust-genotypes --include-homozygous" %(snp_vcf,in_dict['sam_path'], snp_vcf, in_dict['fasta_path'], in_dict['chrom'] ))
            stream.read()
            
            stream=os.popen("bcftools view -e  'GT=\"0\\0\"' %s.phased.preclean.vcf|bgziptabix %s.phased.vcf.gz" %(snp_vcf,snp_vcf))
            stream.read()
            
            stream=os.popen('whatshap --version')
            whatshap_version=stream.read()
            whatshap_version=float(whatshap_version[9:13])
            
            
            if whatshap_version>=0.19:
                stream=os.popen("whatshap haplotag --ignore-read-groups --ignore-linked-read -o %s.phased.bam --reference %s %s.phased.vcf.gz %s --regions %s:%d:%d --ignore-read-groups --ignore-linked-read" %(snp_vcf,in_dict['fasta_path'], snp_vcf, in_dict['sam_path'], args.chrom,start,end))
                stream.read()
            
            else: 
                stream=os.popen('samtools view -b %s %s:%d-%d -o %s.bam --output-fmt BAM; samtools index %s.bam' %(in_dict['sam_path'], args.chrom, start, end, snp_vcf,snp_vcf))
                stream=os.popen("whatshap haplotag --ignore-read-groups --ignore-linked-read -o %s.phased.bam --reference %s %s.phased.vcf.gz %s.bam --ignore-read-groups --ignore-linked-read" %(snp_vcf,in_dict['fasta_path'], snp_vcf, snp_vcf))
                stream.read()
                
            stream=os.popen('samtools index %s.phased.bam' %snp_vcf )
    
    if args.mode in ['indels','both']:
        stream=os.popen('samtools faidx %s %s>%s/%s.fa' %(args.ref,args.chrom,args.vcf,args.chrom))
        stream.read()

        steam=os.popen('rtg format -f fasta %s/%s.fa -o %s/ref.sdf' % (args.chrom,args.vcf,args.vcf))
        stream.read()
        sam_path= '%s.phased.bam' %snp_vcf if args.mode=='both' else args.bam
        
        in_dict={'chrom':args.chrom, 'start':start, 'end':end, 'sam_path':sam_path, 'fasta_path':args.ref, \
             'mincov':args.mincov,  'maxcov':args.maxcov, 'min_allele_freq':args.min_allele_freq, 'min_nbr_sites':args.min_nbr_sites, \
             'threshold':threshold, 'model':args.model, 'cpu':args.cpu,  'vcf_path':args.vcf,'prefix':args.prefix,'sample':args.sample, \
                'del_t':args.del_threshold,'ins_t':args.ins_threshold}
        
        indel_vcf=indelCaller.test_model(in_dict)
        
        if args.mode=='both':
            
            stream=os.popen('rtg vcfdecompose -i %s.vcf.gz --break-mnps --break-indels -o %s.decomposed.vcf.gz -t %s/ref.sdf' %(indel_vcf,indel_vcf,args.vcf))
            stream.read()
            final_path=os.path.join(args.vcf,'%s.final.vcf.gz' %args.prefix)
            stream=os.popen('bcftools concat %s.vcf.gz %s.decomposed.vcf.gz -a -d all |bgziptabix %s' %(snp_vcf, indel_vcf, final_path))
            stream.read()
            
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
