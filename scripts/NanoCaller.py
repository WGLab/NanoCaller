from warnings import simplefilter 
simplefilter(action='ignore', category=FutureWarning)

import time, argparse, os, shutil, sys, pysam, datetime
import multiprocessing as mp
from intervaltree import Interval, IntervalTree
from subprocess import PIPE, Popen
from utils import *

def run(args):
    import snpCaller, indelCaller
        
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
                print('%s: contig %s not found in reference.' %(str(datetime.datetime.now()), args.chrom), flush=True)
                return 
                
        except FileNotFoundError:
            print('%s: Index file .fai required for reference genome file' %(str(datetime.datetime.now())), flush=True)
            return
            
    else:
        end=args.end
    
    if not args.start:
        start=1
    else:
        start=args.start
    
    threshold=[float(args.neighbor_threshold.split(',')[0]), float(args.neighbor_threshold.split(',')[1])]
    
    dirname = os.path.dirname(__file__)
    
    if args.exclude_bed in ['hg38', 'hg19', 'mm10', 'mm39']:
        args.exclude_bed=os.path.join(dirname, 'release_data/bed_files/%s_centro_telo.bed.gz' %args.exclude_bed)
    
    if args.include_bed:
        tbx = pysam.TabixFile(args.include_bed)
        include_intervals=IntervalTree(Interval(int(row[1]), int(row[2]), "%s" % (row[1])) for row in tbx.fetch(args.chrom, parser=pysam.asBed()))
        
        include_intervals=IntervalTree(include_intervals.overlap(start,end))
        
        if include_intervals:
            start=max(start, min(x[0] for x in include_intervals))
            end=min(end, max(x[1] for x in include_intervals))
        
        else:
            print('%s: No overlap between include_bed file and start/end coordinates' %(str(datetime.datetime.now())),flush=True)
            return
        
    in_dict={'chrom':args.chrom, 'start':start, 'end':end, 'sam_path':args.bam, 'fasta_path':args.ref, \
             'mincov':args.mincov,  'maxcov':args.maxcov, 'min_allele_freq':args.min_allele_freq, 'min_nbr_sites':args.min_nbr_sites, \
             'threshold':threshold, 'model':args.model, 'cpu':args.cpu,  'vcf_path':args.vcf,'prefix':args.prefix,'sample':args.sample, \
            'seq':args.sequencing, 'supplementary':args.supplementary, 'include_bed':args.include_bed, 'exclude_bed':args.exclude_bed}
        
    snp_vcf=''
    if args.mode in ['snps', 'snps_unphased', 'both']:
        snp_time=time.time()
        snp_vcf=snpCaller.test_model(in_dict, pool)
        print('\n%s: SNP calling completed for contig %s. Time taken= %.4f\n' %(str(datetime.datetime.now()), in_dict['chrom'], time.time()-snp_time),flush=True)

        if snp_vcf and args.mode in ['snps', 'both']:
            enable_whatshap = '--distrust-genotypes --include-homozygous' if args.enable_whatshap else ''
            
            print('\n%s: ------WhatsHap SNP phasing log------\n' %(str(datetime.datetime.now())),flush=True)
            
            run_cmd("whatshap phase %s.vcf.gz %s -o %s.phased.preclean.vcf -r %s --ignore-read-groups --chromosome %s %s" %(snp_vcf,in_dict['sam_path'], snp_vcf, in_dict['fasta_path'], in_dict['chrom'] , enable_whatshap), verbose=True)
            
            
            run_cmd("bcftools view -e  'GT=\"0\\0\"' %s.phased.preclean.vcf|bgziptabix %s.phased.vcf.gz" %(snp_vcf,snp_vcf))
            
            print('\n%s: ------SNP phasing completed------\n' %(str(datetime.datetime.now())),flush=True)


            if args.mode=='both':
                print('\n%s: ------WhatsHap BAM phasing log------\n' %(str(datetime.datetime.now())),flush=True)
                
                run_cmd("whatshap haplotag --ignore-read-groups --ignore-linked-read -o %s.phased.bam --reference %s %s.phased.vcf.gz %s --regions %s:%d:%d --tag-supplementary" %(snp_vcf,in_dict['fasta_path'], snp_vcf, in_dict['sam_path'], args.chrom,start,end), verbose=True)

                run_cmd('samtools index %s.phased.bam' %snp_vcf )
                
                print('\n%s: ------BAM phasing completed-----\n' %(str(datetime.datetime.now())),flush=True)
            
        else:
            return
            
    if args.mode in ['indels','both']:
        
        sam_path= '%s.phased.bam' %snp_vcf if args.mode=='both' else args.bam
        
        in_dict={'chrom':args.chrom, 'start':start, 'end':end, 'sam_path':sam_path, 'fasta_path':args.ref, \
             'mincov':args.mincov,  'maxcov':args.maxcov, 'min_allele_freq':args.min_allele_freq, 'min_nbr_sites':args.min_nbr_sites, \
             'threshold':threshold, 'model':args.model, 'cpu':args.cpu,  'vcf_path':args.vcf,'prefix':args.prefix,'sample':args.sample, 'seq':args.sequencing, \
                'del_t':args.del_threshold,'ins_t':args.ins_threshold,'supplementary':args.supplementary, 'include_bed':args.include_bed\
                , 'exclude_bed':args.exclude_bed}
        ind_time=time.time()
        indel_vcf=indelCaller.test_model(in_dict, pool)
        
        print('%s: Indel calling completed for contig %s. Time taken= %.4f' %(str(datetime.datetime.now()), in_dict['chrom'], time.time()-ind_time),flush=True)
        
        if args.mode=='both':
            print('%s: Post processing' %(str(datetime.datetime.now())),flush=True)
            run_cmd('samtools faidx %s %s>%s/%s.fa' %(args.ref,args.chrom,args.vcf,args.chrom))

            if os.path.exists('%s/ref.sdf' %args.vcf):
                if os.path.isdir('%s/ref.sdf' %args.vcf):
                    shutil.rmtree('%s/ref.sdf' %args.vcf)
                else:
                    os.remove('%s/ref.sdf' %args.vcf)

            run_cmd('rtg RTG_MEM=4G format -f fasta %s/%s.fa -o %s/ref.sdf' % (args.vcf,args.chrom,args.vcf))

            if os.path.exists('%s.decomposed.vcf.gz' %indel_vcf):
                if os.path.isdir('%s.decomposed.vcf.gz' %indel_vcf):
                    shutil.rmtree('%s.decomposed.vcf.gz' %indel_vcf)
                else:
                    os.remove('%s.decomposed.vcf.gz' %indel_vcf)
                
                
            run_cmd('rtg RTG_MEM=4G vcfdecompose -i %s.vcf.gz --break-mnps -o - -t %s/ref.sdf|rtg vcffilter -i - --non-snps-only -o  %s.decomposed.vcf.gz' %(indel_vcf,args.vcf,indel_vcf))

            final_path=os.path.join(args.vcf,'%s.final.vcf.gz' %args.prefix)
            run_cmd('bcftools concat %s.phased.vcf.gz %s.decomposed.vcf.gz -a -d all |bgziptabix %s' %(snp_vcf, indel_vcf, final_path))

    pool.close()
    pool.join()
    
if __name__ == '__main__':
    
    t=time.time()
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('Required arguments')
    parser.add_argument("-mode",  "--mode",  help="Testing mode, options are 'snps', 'snps_unphased', 'indels' and 'both'. 'snps_unphased' mode quits NanoCaller without using WhatsHap for phasing.", type=str, default='both')
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
    
    parser.add_argument("-include_bed",  "--include_bed",  help="Only call variants inside the intervals specified in the bgzipped and tabix indexed BED file. If any other flags are used to specify a region, intersect the region with intervals in the BED file, e.g. if -chom chr1 -start 10000000 -end 20000000 flags are set, call variants inside the intervals specified by the BED file that overlap with chr1:10000000-20000000. Same goes for the case when whole genome variant calling flag is set.", type=str, default=None)
    
    parser.add_argument("-exclude_bed",  "--exclude_bed",  help="Path to bgzipped and tabix indexed BED file containing intervals to ignore  for variant calling. BED files of centromere and telomere regions for the following genomes are included in NanoCaller: hg38, hg19, mm10 and mm39. To use these BED files use one of the following options: ['hg38', 'hg19', 'mm10', 'mm39'].", type=str, default=None)
    
    parser.add_argument("-sample",  "--sample",  help="VCF file sample name", type=str, default='SAMPLE')
    parser.add_argument("-sup",  "--supplementary",  help="Use supplementary reads", default=False, action='store_true')
    parser.add_argument("-mincov",  "--mincov",  help="min coverage", type=int, default=8)
    parser.add_argument("-maxcov",  "--maxcov",  help="max coverage", type=int, default=160)
    
    parser.add_argument("-start",  "--start",  help="start, default is 1", type=int)
    parser.add_argument("-end",  "--end",  help="end, default is the end of contig", type=int)
    parser.add_argument("-nbr_t",  "--neighbor_threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended and for PacBio reads '0.3,0.7' is recommended", type=str, default='0.4,0.6')
    parser.add_argument("-ins_t", "--ins_threshold", help="Insertion Threshold",type=float,default=0.4)
    parser.add_argument("-del_t", "--del_threshold", help="Deletion Threshold",type=float,default=0.6)
        
    parser.add_argument("-enable_whatshap",  "--enable_whatshap",  help="Allow WhatsHap to change SNP genotypes when phasing using --distrust-genotypes and --include-homozygous flags (this is not the same as regenotyping), considerably increasing the time needed for phasing. It has a negligible effect on SNP calling accuracy for Nanopore reads, but may make a small improvement for PacBio reads. By default WhatsHap will only phase SNP calls produced by NanoCaller, but not change their genotypes.",  default=False, action='store_true')
    
    parser.add_argument('-wgs_print_commands','--wgs_print_commands', help='If set, print the commands to run NanoCaller on all contigs in a file named "wg_commands". By default, run the NanoCaller on each contig in a sequence.', default=False, action='store_true')
    
    parser.add_argument('-wgs_contigs_type','--wgs_contigs_type', \
                        help="""Options are "with_chr", "without_chr" and "all",\
                        or a space/whitespace separated list of contigs in quotation\
                        marks e.g. "chr3 chr6 chr22" . "with_chr" option will assume \
                        human genome and run NanoCaller on chr1-22, "without_chr" will \
                        run on chromosomes 1-22 if the BAM and reference genome files \
                        use chromosome names without "chr". "all" option will run \
                        NanoCaller on each contig present in reference genome FASTA file.""", \
                        type=str, default='with_chr')
    
    
    
    args = parser.parse_args()
    
    if not args.vcf:
        args.vcf=os.getcwd()
    
    try:  
        os.mkdir(args.vcf)  
    except OSError as error:  
        pass
    
    with open(os.path.join(args.vcf,'args'),'w') as file:
        file.write(str(args))
        
    print('\n%s: Starting NanoCaller.\n\nRunning arguments are saved in the following file: %s\n' %(str(datetime.datetime.now()), os.path.join(args.vcf,'args')), flush=True)
    
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
                print('%s: index file .fai required for reference genome file' %str(datetime.datetime.now()), flush=True)
                sys.exit()
            
        else:
            chrom_list= args.wgs_contigs_type.split()
        
        if args.include_bed:
            stream=run_cmd('zcat %s|cut -f 1|uniq' %args.include_bed, output=True)
            bed_chroms=stream.split()
            
            chrom_list=[chrom for chrom in chrom_list if chrom in bed_chroms]
            
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
            print('%s: Commands for running NanoCaller on contigs in whole genome are saved in the file %s' %(str(datetime.datetime.now()), os.path.join(args.vcf,'wg_commands')))
        
        else:
            vcf_path=args.vcf
            for chrom in chrom_list:
                ctime=time.time()
                args.chrom=chrom
                args.vcf=os.path.join(vcf_path, args.chrom)
                run(args)
                print('%s: Variant calling on %s completed. Time taken=%.4fs.' %(str(datetime.datetime.now()), chrom, time.time()-ctime))
        
    elapsed=time.time()-t
    print ('\n%s: Total Time Elapsed: %.2f seconds' %(str(datetime.datetime.now()), elapsed))
