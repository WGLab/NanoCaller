import time, argparse, os, shutil, sys
import multiprocessing as mp

def run(args):
    pool = mp.Pool(processes=args.cpu)

    if not args.vcf:
        args.vcf=os.getcwd()
    
    try:  
        os.mkdir(args.vcf)  
    except OSError as error:  
        pass 
    
    
    with open(os.path.join(args.vcf,'args'),'w') as file:
        file.write(str(args))
        
    
    
    end=None
    if not args.end:
        try:
            with open(args.ref+'.fai','r') as file:
                for line in file:
                    if line.split('\t')[0]==args.chrom:
                        
                        end=int(line.split('\t')[1])
            
            if end==None:
                print('contig %s not found in reference' %args.chrom, flush=True)
                return 
                
        except FileNotFoundError:
            print('index file .fai required for reference genome file', flush=True)
            return
            
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
            'seq':args.sequencing, 'supplementary':args.supplementary}
        
    snp_vcf=''
    if args.mode in ['snps','both']:
        snp_time=time.time()
        snp_vcf=snpCaller.test_model(in_dict, pool)
        print('SNP calling completed for contig %s. Time taken= %.4f' %(in_dict['chrom'], time.time()-snp_time),flush=True)

        if snp_vcf:
            disable_whatshap = ' --distrust-genotypes --include-homozygous' if args.disable_whatshap else ''
            stream=os.popen("whatshap phase %s.vcf.gz %s -o %s.phased.preclean.vcf -r %s --ignore-read-groups --chromosome %s %s" %(snp_vcf,in_dict['sam_path'], snp_vcf, in_dict['fasta_path'], in_dict['chrom'] , disable_whatshap))
            stream.read()
            
            stream=os.popen("bcftools view -e  'GT=\"0\\0\"' %s.phased.preclean.vcf|bgziptabix %s.phased.vcf.gz" %(snp_vcf,snp_vcf))
            stream.read()


            if args.mode=='both':
                stream=os.popen("whatshap haplotag --ignore-read-groups --ignore-linked-read -o %s.phased.bam --reference %s %s.phased.vcf.gz %s --regions %s:%d:%d --tag-supplementary" %(snp_vcf,in_dict['fasta_path'], snp_vcf, in_dict['sam_path'], args.chrom,start,end))
                stream.read()


                stream=os.popen('samtools index %s.phased.bam' %snp_vcf )
                stream.read()
            
        else:
            return
            
    if args.mode in ['indels','both']:
        
        sam_path= '%s.phased.bam' %snp_vcf if args.mode=='both' else args.bam
        
        in_dict={'chrom':args.chrom, 'start':start, 'end':end, 'sam_path':sam_path, 'fasta_path':args.ref, \
             'mincov':args.mincov,  'maxcov':args.maxcov, 'min_allele_freq':args.min_allele_freq, 'min_nbr_sites':args.min_nbr_sites, \
             'threshold':threshold, 'model':args.model, 'cpu':args.cpu,  'vcf_path':args.vcf,'prefix':args.prefix,'sample':args.sample, 'seq':args.sequencing, \
                'del_t':args.del_threshold,'ins_t':args.ins_threshold,'supplementary':args.supplementary}
        ind_time=time.time()
        indel_vcf=indelCaller.test_model(in_dict, pool)
        print('Indel calling completed for contig %s. Time taken= %.4f' %(in_dict['chrom'], time.time()-ind_time),flush=True)
        
        if args.mode=='both':
            print('Post processing',flush=True)
            stream=os.popen('samtools faidx %s %s>%s/%s.fa' %(args.ref,args.chrom,args.vcf,args.chrom))
            stream.read()

            if os.path.exists('%s/ref.sdf' %args.vcf):
                if os.path.isdir('%s/ref.sdf' %args.vcf):
                    shutil.rmtree('%s/ref.sdf' %args.vcf)
                else:
                    os.remove('%s/ref.sdf' %args.vcf)

            stream=os.popen('rtg RTG_MEM=4G format -f fasta %s/%s.fa -o %s/ref.sdf' % (args.vcf,args.chrom,args.vcf))
            stream.read()

            if os.path.exists('%s.decomposed.vcf.gz' %indel_vcf):
                if os.path.isdir('%s.decomposed.vcf.gz' %indel_vcf):
                    shutil.rmtree('%s.decomposed.vcf.gz' %indel_vcf)
                else:
                    os.remove('%s.decomposed.vcf.gz' %indel_vcf)
                
                
            stream=os.popen('rtg RTG_MEM=4G vcfdecompose -i %s.vcf.gz --break-mnps --break-indels -o %s.decomposed.vcf.gz -t %s/ref.sdf' %(indel_vcf,indel_vcf,args.vcf))
            stream.read()

            final_path=os.path.join(args.vcf,'%s.final.vcf.gz' %args.prefix)
            stream=os.popen('bcftools concat %s.phased.vcf.gz %s.decomposed.vcf.gz -a -d all |bgziptabix %s' %(snp_vcf, indel_vcf, final_path))
            stream.read()

    pool.close()
    pool.join()
    
if __name__ == '__main__':
    t=time.time()
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('Required arguments')
    parser.add_argument("-mode",  "--mode",  help="Testing mode, options are 'snps', 'indels' and 'both'", type=str, default='both')
    parser.add_argument("-seq",  "--sequencing",  help="Sequencing type, options are 'ont' and 'pacbio'", type=str, default='ont')
    parser.add_argument("-model",  "--model",  help="NanoCaller SNP model to be used, options are 'NanoCaller1' (trained on HG001 Nanopore reads), 'NanoCaller2' (trained on HG002 Nanopore reads) and 'NanoCaller3' (trained on HG003 PacBio reads) ", default='NanoCaller1')
    parser.add_argument("-vcf",  "--vcf",  help="VCF output path, default is current working directory", type=str)
    
    parser.add_argument("-chrom",  "--chrom",  help="Chromosome")
    parser.add_argument("-cpu",  "--cpu",  help="CPUs", type=int, default=1)
    parser.add_argument("-min_allele_freq",  "--min_allele_freq",  help="minimum alternative allele frequency", type=float,  default=0.15)
    parser.add_argument("-min_nbr_sites",  "--min_nbr_sites",  help="minimum number of nbr sites", type=int,  default =1)
    
    requiredNamed.add_argument("-bam",  "--bam",  help="Bam file, should be phased if 'indel' mode is selected", required=True)
    requiredNamed.add_argument("-ref",  "--ref",  help="reference genome file with .fai index", required=True)
    
    requiredNamed.add_argument("-prefix",  "--prefix",  help="VCF file prefix", type=str, required=True)
    parser.add_argument("-sample",  "--sample",  help="VCF file sample name", type=str, default='SAMPLE')
    parser.add_argument("-sup",  "--supplementary",  help="Use supplementary reads", default=False, action='store_true')
    parser.add_argument("-mincov",  "--mincov",  help="min coverage", type=int, default=8)
    parser.add_argument("-maxcov",  "--maxcov",  help="max coverage", type=int, default=160)
    
    parser.add_argument("-start",  "--start",  help="start, default is 1", type=int)
    parser.add_argument("-end",  "--end",  help="end, default is the end of contig", type=int)
    parser.add_argument("-nbr_t",  "--neighbor_threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended and for PacBio reads '0.3,0.7' is recommended", type=str, default='0.4,0.6')
    parser.add_argument("-ins_t", "--ins_threshold", help="Insertion Threshold",type=float,default=0.4)
    parser.add_argument("-del_t", "--del_threshold", help="Deletion Threshold",type=float,default=0.6)
    
    parser.add_argument("-disable_whatshap",  "--disable_whatshap",  help="Allow WhatsHap to change SNP genotypes when phasing",  default=False, action='store_true')
    
    parser.add_argument('-wgs_print_commands','--wgs_print_commands', help='If set, print the commands to run NanoCaller on all contigs in a file named "wg_commands". By default, run the NanoCaller on each contig in a sequence.', default=False, action='store_true')
    
    parser.add_argument('-wgs_contigs_type','--wgs_contigs_type', \
                        help='Options are "with_chr", "without_chr" and "all",\ 
                        or a space/whitespace separated list of contigs in quotation\
                        marks e.g. "chr3 chr6 chr22" . "with_chr" option will assume \
                        human genome and run NanoCaller on chr1-22, "without_chr" will \
                        run on chromosomes 1-22 if the BAM and reference genome files \
                        use chromosome names without "chr". "all" option will run \
                        NanoCaller on each contig present in reference genome FASTA file.', \
                        type=str, default='with_chr')
    
    import snpCaller, indelCaller
    
    args = parser.parse_args()
    
    if not args.vcf:
        args.vcf=os.getcwd()
    
    try:  
        os.mkdir(args.vcf)  
    except OSError as error:  
        pass
    
    if args.chrom:
        run(args)
    
    else:
        if args.wgs_contigs_type=='with_chr':
            chrom_list=['chr%d' %d for d in range(1,23)]
        
        elif args.wgs_contigs_type == 'without_chr':
            chrom_list=['%d' %d for d in range(1,23)]
        
        elif args.wgs_contigs_type == 'all':
            chrom_list=[]
            
            try:
                with open(args.ref+'.fai','r') as file:
                    for line in file:
                        chrom_list.append(line.split('\t')[0])
                
            except FileNotFoundError:
                print('index file .fai required for reference genome file', flush=True)
                sys.exit()
            
        else:
            chrom_list= args.wgs_contigs_type.split()
        
        if args.wgs_print_commands:
            args_dict=vars(args)
            
            with open(os.path.join(args.vcf,'wg_commands'),'w') as wg_commands:
                for chrom in chrom_list:
                    cmd=''
                    vcf=os.path.join(args.vcf, chrom)
                    for x in args_dict:
                        if x not in ['chrom','wgs_print_commands','wgs_contigs_type','start','end','vcf']:
                            cmd+= '--%s %s ' %(x, args_dict[x])
                            
                    dirname = os.path.dirname(__file__)
                    wg_commands.write('python %s/NanoCaller.py -chrom %s %s -vcf %s \n' %(dirname, chrom, cmd, vcf))
            print('Commands for running NanoCaller on contigs in whole genome are saved in the file %s' %os.path.join(args.vcf,'wg_commands'))
        
        else:
            vcf_path=args.vcf
            for chrom in chrom_list:
                ctime=time.time()
                args.chrom=chrom
                args.vcf=os.path.join(vcf_path, args.chrom)
                run(args)
                print('Variant calling on %s completed. Time taken=%.4fs.' %(chrom,time.time()-ctime))
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
