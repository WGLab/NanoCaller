from warnings import simplefilter 
simplefilter(action='ignore', category=FutureWarning)

import time, argparse, os, shutil, sys, pysam, datetime
import multiprocessing as mp
from intervaltree import Interval, IntervalTree
from subprocess import PIPE, Popen
from utils import *
    
if __name__ == '__main__':
    
    t=time.time()
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('Required arguments')
    parser.add_argument("-mode",  "--mode",  help="Testing mode, options are 'snps', 'snps_unphased', 'indels' and 'both'. 'snps_unphased' mode quits NanoCaller without using WhatsHap for phasing.", type=str, default='both')
    parser.add_argument("-seq",  "--sequencing",  help="Sequencing type, options are 'ont' and 'pacbio'", type=str, default='ont')
    parser.add_argument("-model",  "--model",  help="NanoCaller SNP model to be used, options are 'NanoCaller1' (trained on HG001 Nanopore reads), 'NanoCaller2' (trained on HG002 Nanopore reads) and 'NanoCaller3' (trained on HG003 PacBio reads) ", default='NanoCaller1')
    parser.add_argument("-o",  "--output",  help="VCF output path, default is current working directory", type=str)
    
    parser.add_argument("-chrom",  "--chrom",  help='A space/whitespace separated list of contigs in quotation marks e.g. "chr3 chr6 chr22" .')
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
    parser.add_argument("-nbr_t",  "--neighbor_threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended, for PacBio CCS anc CLR reads '0.3,0.7' and '0.3,0.6' are recommended respectively", type=str, default='0.4,0.6')
    parser.add_argument("-ins_t", "--ins_threshold", help="Insertion Threshold",type=float,default=0.4)
    parser.add_argument("-del_t", "--del_threshold", help="Deletion Threshold",type=float,default=0.6)
        
    parser.add_argument("-enable_whatshap",  "--enable_whatshap",  help="Allow WhatsHap to change SNP genotypes when phasing using --distrust-genotypes and --include-homozygous flags (this is not the same as regenotyping), considerably increasing the time needed for phasing. It has a negligible effect on SNP calling accuracy for Nanopore reads, but may make a small improvement for PacBio reads. By default WhatsHap will only phase SNP calls produced by NanoCaller, but not change their genotypes.",  default=False, action='store_true')
    
    parser.add_argument('-keep_bam','--keep_bam', help='Keep phased bam files.', default=False, action='store_true')
    
    
    parser.add_argument('-wgs_contigs_type','--wgs_contigs_type', \
                        help="""Options are "with_chr", "without_chr" and "all",\
                        "with_chr" option will assume \
                        human genome and run NanoCaller on chr1-22, "without_chr" will \
                        run on chromosomes 1-22 if the BAM and reference genome files \
                        use chromosome names without "chr". "all" option will run \
                        NanoCaller on each contig present in reference genome FASTA file.""", \
                        type=str, default='with_chr')
    
    
    
    args = parser.parse_args()
    
    if not args.output:
        args.output=os.getcwd()
    
    os.makedirs(args.output, exist_ok=True)  
    
    with open(os.path.join(args.output,'args'),'w') as file:
        file.write(str(args))
        
    print('\n%s: Starting NanoCaller.\n\nRunning arguments are saved in the following file: %s\n' %(str(datetime.datetime.now()), os.path.join(args.output,'args')), flush=True)
    
    if args.chrom:
        chrom_list= args.chrom.split()
    
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

    if args.include_bed:
        stream=run_cmd('zcat %s|cut -f 1|uniq' %args.include_bed, output=True)
        bed_chroms=stream.split()

        chrom_list=[chrom for chrom in chrom_list if chrom in bed_chroms]

    args_dict=vars(args)
    chrom_lengths={}
    with open(args.ref+'.fai','r') as file:
        for line in file:
            chrom_lengths[line.split('\t')[0]]=int(line.split('\t')[1])

    with open(os.path.join(args.output,'wg_commands'),'w') as wg_commands:
        job_counter=0
        for chrom in chrom_list:
            cmd=''

            for x in args_dict:
                if x in ['chrom','wgs_contigs_type','start','end','output','cpu','prefix'] or args_dict[x] is None:
                    pass

                elif x in ['supplementary', 'enable_whatshap','keep_bam']:
                    if args_dict[x]==True:
                        cmd+=' --%s ' %x

                else:
                    cmd+= '--%s %s ' %(x, args_dict[x])

            dirname = os.path.dirname(__file__)

            try:
               chr_end=chrom_lengths[chrom]
            
            except KeyError:
               print('Contig %s not found in reference' %chrom,flush=True)
               continue
            for mbase in range(1,chr_end,10000000):
                out_path=os.path.join(args.output, 'intermediate_files', '%s_%d_%d' %(chrom, mbase, min(chr_end,mbase+10000000-1)))
                wg_commands.write('python %s -chrom %s %s -cpu 1 --output %s -start %d -end %d -prefix %s_%d_%d\n' %(os.path.join(dirname,'NanoCaller.py'), chrom, cmd, out_path ,mbase, min(chr_end,mbase+10000000-1),chrom, mbase, min(chr_end,mbase+10000000-1)))
                
                job_counter+=1
                
    print('%s: Commands for running NanoCaller on contigs in whole genome are saved in the file %s' %(str(datetime.datetime.now()), os.path.join(args.output,'wg_commands')))
    
    print('Running %d jobs using %d workers in parallel.\n' %(job_counter, args.cpu))

    run_cmd('cat %s|parallel -j %d' %(os.path.join(args.output,'wg_commands'), args.cpu))
    
    out_path=os.path.join(args.output, 'intermediate_files')
                 
    if args.mode in ['snps_unphased','snps','both']:
        run_cmd('ls -1 %s/*/*snps.vcf.gz|bcftools concat -f -|bcftools sort|bgziptabix %s/%s.snps.vcf.gz' %(out_path, args.output, args.prefix))
        
        if args.mode!='snps_unphased':
            run_cmd('ls -1 %s/*/*snps.phased.vcf.gz|bcftools concat -f -|bcftools sort|bgziptabix %s/%s.snps.phased.vcf.gz' %(out_path, args.output, args.prefix))
    
    if args.mode in ['indels','both']:
        run_cmd('ls -1 %s/*/*indels.vcf.gz|bcftools concat -f -|bcftools sort|bgziptabix %s/%s.indels.vcf.gz' %(out_path, args.output, args.prefix))
    
    if args.mode=='both':
        run_cmd('ls -1 %s/*/*final.vcf.gz|bcftools concat -f -|bcftools sort|bgziptabix %s/%s.final.vcf.gz' %(out_path, args.output, args.prefix))
    
    
    elapsed=time.time()-t
    print ('\n%s: Total Time Elapsed: %.2f seconds' %(str(datetime.datetime.now()), elapsed))
