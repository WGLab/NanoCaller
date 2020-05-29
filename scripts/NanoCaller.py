import time, argparse, os, snpCaller, indelCaller

if __name__ == '__main__':
    t=time.time()
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-mode",  "--mode",  help="Testing mode, options are 'snps', 'indels' and 'both'",type=str,default='both')
    parser.add_argument("-model",  "--model",  help="Model")
    parser.add_argument("-vcf",  "--vcf",  help="VCF output path",type=str,default='')
    
    parser.add_argument("-chrom",  "--chrom",  help="Chromosome")
    parser.add_argument("-cpu",  "--cpu",  help="CPUs", type=int,default=1)
    parser.add_argument("-min_allele_freq",  "--min_allele_freq",  help="minimum alternative allele frequency", type=float,  default=0.15)
    parser.add_argument("-min_nbr_sites",  "--min_nbr_sites",  help="minimum number of nbr sites", type=int,  default =1)
    
    parser.add_argument("-bam",  "--bam",  help="Bam file")
    parser.add_argument("-ref",  "--ref",  help="reference genome file with .fai index")
    
    parser.add_argument("-prefix",  "--prefix",  help="VCF file prefix",type=str,default='')
    parser.add_argument("-sample",  "--sample",  help="VCF file sample name",type=str,default='SAMPLE')
    parser.add_argument("-bed",  "--bed",  help="BED file")
    
    parser.add_argument("-mincov",  "--mincov",  help="min coverage", type=int, default=8)
    parser.add_argument("-maxcov",  "--maxcov",  help="max coverage", type=int, default=160)
    
    parser.add_argument("-start",  "--start",  help="start, default is 1", type=int)
    parser.add_argument("-end",  "--end",  help="end, default is the end of contig", type=int)
    parser.add_argument("-threshold",  "--threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, default=0.3,0.7", type=str, default='0.3,0.7')
    parser.add_argument("-ins_t", "--ins_threshold", help="Insertion Threshold",type=float,default=0.4)
    parser.add_argument("-del_t", "--del_threshold", help="Deletion Threshold",type=float,default=0.6)
    parser.add_argument("-remove_homopolymer", "--remove_homopolymer", help="Remove Homopolymer regions from indel calling", type=bool, default=True)
    
    args = parser.parse_args()
    
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
    
    threshold=[float(args.threshold.split(',')[0]), float(args.threshold.split(',')[1])]
    
    in_dict={'chrom':args.chrom, 'start':start, 'end':end, 'sam_path':args.bam, 'fasta_path':args.ref, 'bed':args.bed, \
             'mincov':args.mincov,  'maxcov':args.maxcov, 'min_allele_freq':args.min_allele_freq, 'min_nbr_sites':args.min_nbr_sites, \
             'threshold':threshold, 'model':args.model, 'cpu':args.cpu,  'vcf_path':args.vcf,'prefix':args.prefix,'sample':args.sample}
        
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
        sam_path= '%s.phased.bam' %snp_vcf if args.mode=='both' else args.bam
        
        in_dict={'chrom':args.chrom, 'start':start, 'end':end, 'sam_path':sam_path, 'fasta_path':args.ref, 'bed':args.bed, \
             'mincov':args.mincov,  'maxcov':args.maxcov, 'min_allele_freq':args.min_allele_freq, 'min_nbr_sites':args.min_nbr_sites, \
             'threshold':threshold, 'model':args.model, 'cpu':args.cpu,  'vcf_path':args.vcf,'prefix':args.prefix,'sample':args.sample, \
                'del_t':args.del_threshold,'ins_t':args.ins_threshold,'remove_homopolymer':args.remove_homopolymer}
        
        indelCaller.test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)