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

    if not args.output:
        args.output=os.getcwd()
    
    os.makedirs(args.output, exist_ok=True)
    
    end=None
    if not args.end:
        try:
            with open(args.ref+'.fai','r') as file:
                for line in file:
                    if line.split('\t')[0]==args.chrom:
                        
                        end=int(line.split('\t')[1])
            
            if end==None:
                print('%s: contig %s not found in reference.' %(str(datetime.datetime.now()), args.chrom), flush=True)
                sys.exit(2) 
                
        except FileNotFoundError:
            print('%s: Index file .fai required for reference genome file' %(str(datetime.datetime.now())), flush=True)
            sys.exit(2)
            
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
             'threshold':threshold, 'snp_model':args.snp_model, 'cpu':args.cpu,  'vcf_path':args.output,'prefix':args.prefix,'sample':args.sample, \
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


            if args.mode=='both' or args.phase_bam:
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
             'threshold':threshold, 'snp_model':args.snp_model,'indel_model':args.indel_model, 'cpu':args.cpu,  'vcf_path':args.output,'prefix':args.prefix,'sample':args.sample, 'seq':args.sequencing, \
                'del_t':args.del_threshold,'ins_t':args.ins_threshold,'supplementary':args.supplementary, 'include_bed':args.include_bed\
                , 'exclude_bed':args.exclude_bed,'win_size':args.win_size,'small_win_size':args.small_win_size}
        ind_time=time.time()
        indel_vcf=indelCaller.test_model(in_dict, pool)
        
        print('%s: Post processing' %(str(datetime.datetime.now())),flush=True)
        
        run_cmd('samtools faidx %s %s > %s/%s.fa' %(args.ref,args.chrom,args.output,args.chrom))

        remove_path('%s/ref.sdf' %args.output)

        run_cmd('rtg RTG_MEM=4G format -f fasta %s/%s.fa -o %s/ref.sdf' % (args.output,args.chrom,args.output))

        remove_path('%s.vcf.gz' %indel_vcf)    

        run_cmd('rtg RTG_MEM=4G vcfdecompose -i %s.raw.vcf.gz --break-mnps -o - -t %s/ref.sdf|rtg RTG_MEM=4G vcffilter -i - --non-snps-only -o  %s.vcf.gz' %(indel_vcf,args.output,indel_vcf))
            
        print('%s: Indel calling completed for contig %s. Time taken= %.4f' %(str(datetime.datetime.now()), in_dict['chrom'], time.time()-ind_time),flush=True)
        
        if args.mode=='both':
            
            if not args.keep_bam and not args.phase_bam:
                os.remove('%s.phased.bam' %snp_vcf)
                
            final_path=os.path.join(args.output,'%s.final.vcf.gz' %args.prefix)
            run_cmd('bcftools concat %s.phased.vcf.gz %s.vcf.gz -a -d all |bgziptabix %s' %(snp_vcf, indel_vcf, final_path))

    pool.close()
    pool.join()
    
if __name__ == '__main__':
    
    t=time.time()
    
    
    preset_dict={'ont':{'sequencing':'ont', 'snp_model':'ONT-HG002', 'indel_model':'ONT-HG002', 'neighbor_threshold':'0.4,0.6', 'ins_threshold':0.4,'del_threshold':0.6, 'enable_whatshap':False},
                 
                'ul_ont': {'sequencing':'ul_ont', 'snp_model':'ONT-HG002', 'indel_model':'ONT-HG002', 'neighbor_threshold':'0.4,0.6', 'ins_threshold':0.4,'del_threshold':0.6, 'enable_whatshap':False},
                 
                'ul_ont_extreme':{'sequencing':'ul_ont_extreme', 'snp_model':'ONT-HG002', 'indel_model':'ONT-HG002', 'neighbor_threshold':'0.4,0.6', 'ins_threshold':0.4,'del_threshold':0.6, 'enable_whatshap':False},
                 
                'ccs':{'sequencing':'pacbio', 'snp_model': 'CCS-HG002', 'indel_model':'CCS-HG002', 'neighbor_threshold':'0.3,0.7', 'ins_threshold':0.4,'del_threshold':0.4, 'enable_whatshap':True},
                 
                'clr':{'sequencing':'pacbio', 'snp_model':'CLR-HG002', 'indel_model':'ONT-HG002', 'neighbor_threshold':'0.3,0.6', 'ins_threshold':0.6,'del_threshold':0.6, 'win_size':10, 'small_win_size':2, 'enable_whatshap':True}
                }
    
    flag_dict={"seq":"sequencing", "p":"preset", "o":"output", "sup":"supplementary","nbr_t":"neighbor_threshold","ins_t":"ins_threshold", "del_t":"del_threshold"}
    flag_map=lambda x: flag_dict[x] if x in flag_dict else x
    
    
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    requiredNamed = parser.add_argument_group('Required Arguments')
    preset_group=parser.add_argument_group("Preset")
    config_group = parser.add_argument_group('Configurations')
    region_group=parser.add_argument_group("Variant Calling Regions")
    
    snp_group=parser.add_argument_group("SNP Calling")
    indel_group=parser.add_argument_group("Indel Calling")
    out_group=parser.add_argument_group("Output Options")
    phase_group=parser.add_argument_group("Phasing")
    
    config_group.add_argument("-mode",  "--mode",  help="NanoCaller mode to run, options are 'snps', 'snps_unphased', 'indels' and 'both'. 'snps_unphased' mode quits NanoCaller without using WhatsHap for phasing.", type=str, default='both')
    config_group.add_argument("-seq",  "--sequencing",  help="Sequencing type, options are 'ont', 'ul_ont', 'ul_ont_extreme', and 'pacbio'.  'ont' works well for any type of ONT sequencing datasets. However, use 'ul_ont' if you have several ultra-long ONT reads up to 100kbp long, and 'ul_ont_extreme' if you have several ultra-long ONT reads up to 300kbp long. For PacBio CCS (HiFi) and CLR reads, use 'pacbio'.", type=str, default='ont')

    config_group.add_argument("-cpu",  "--cpu",  help="Number of CPUs to use", type=int, default=1)
    
    config_group.add_argument("-mincov",  "--mincov",  help="Minimum coverage to call a variant", type=int, default=8)
    config_group.add_argument("-maxcov",  "--maxcov",  help="Maximum coverage of reads to use. If sequencing depth at a candidate site exceeds maxcov then reads are downsampled.", type=int, default=160)
    
    #output options
    out_group.add_argument('-keep_bam','--keep_bam', help='Keep phased bam files.', default=False, action='store_true')
    out_group.add_argument("-o",  "--output",  help="VCF output path, default is current working directory", type=str)    
    out_group.add_argument("-prefix",  "--prefix",  help="VCF file prefix", type=str, default='variant_calls')
    out_group.add_argument("-sample",  "--sample",  help="VCF file sample name", type=str, default='SAMPLE')

    #region
    region_group.add_argument("-include_bed",  "--include_bed",  help="Only call variants inside the intervals specified in the bgzipped and tabix indexed BED file. If any other flags are used to specify a region, intersect the region with intervals in the BED file, e.g. if -chom chr1 -start 10000000 -end 20000000 flags are set, call variants inside the intervals specified by the BED file that overlap with chr1:10000000-20000000. Same goes for the case when whole genome variant calling flag is set.", type=str, default=None)
    region_group.add_argument("-exclude_bed",  "--exclude_bed",  help="Path to bgzipped and tabix indexed BED file containing intervals to ignore  for variant calling. BED files of centromere and telomere regions for the following genomes are included in NanoCaller: hg38, hg19, mm10 and mm39. To use these BED files use one of the following options: ['hg38', 'hg19', 'mm10', 'mm39'].", type=str, default=None)
    region_group.add_argument("-start",  "--start",  help="start, default is 1", type=int)
    region_group.add_argument("-end",  "--end",  help="end, default is the end of contig", type=int)
    
    #preset
    preset_group.add_argument("-p",  "--preset",  help="Apply recommended preset values for SNP and Indel calling parameters, options are 'ont', 'ul_ont', 'ul_ont_extreme', 'ccs' and 'clr'. 'ont' works well for any type of ONT sequencing datasets. However, use 'ul_ont' if you have several ultra-long ONT reads up to 100kbp long, and 'ul_ont_extreme' if you have several ultra-long ONT reads up to 300kbp long. For PacBio CCS (HiFi) and CLR reads, use 'ccs'and 'clr' respectively. Presets are described in detail here: github.com/WGLab/NanoCaller/blob/master/docs/Usage.md#preset-options.", type=str)
    
    #required
    requiredNamed.add_argument("-bam",  "--bam",  help="Bam file, should be phased if 'indel' mode is selected", required=True)
    requiredNamed.add_argument("-ref",  "--ref",  help="Reference genome file with .fai index", required=True)
    requiredNamed.add_argument("-chrom",  "--chrom",  help="Chromosome",required=True)
    
    #snp
    snp_group.add_argument("-snp_model",  "--snp_model",  help="NanoCaller SNP model to be used", default='ONT-HG002')
    snp_group.add_argument("-min_allele_freq",  "--min_allele_freq",  help="minimum alternative allele frequency", type=float,  default=0.15)
    snp_group.add_argument("-min_nbr_sites",  "--min_nbr_sites",  help="minimum number of nbr sites", type=int,  default =1)
    snp_group.add_argument("-nbr_t",  "--neighbor_threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended, for PacBio CCS anc CLR reads '0.3,0.7' and '0.3,0.6' are recommended respectively", type=str, default='0.4,0.6')
    snp_group.add_argument("-sup",  "--supplementary",  help="Use supplementary reads", default=False, action='store_true')
    
    
    #indel
    indel_group.add_argument("-indel_model",  "--indel_model",  help="NanoCaller indel model to be used", default='ONT-HG002')
    indel_group.add_argument("-ins_t", "--ins_threshold", help="Insertion Threshold",type=float,default=0.4)
    indel_group.add_argument("-del_t", "--del_threshold", help="Deletion Threshold",type=float,default=0.6)
    indel_group.add_argument("-win_size",  "--win_size",  help="Size of the sliding window in which the number of indels is counted to determine indel candidate site. Only indels longer than 2bp are counted in this window. Larger window size can increase recall, but use a maximum of 50 only", type=int, default=40)
    indel_group.add_argument("-small_win_size",  "--small_win_size",  help="Size of the sliding window in which indel frequency is determined for small indels", type=int, default=4)
    
    
    #phasing
    phase_group.add_argument('-phase_bam','--phase_bam', help='Phase bam files if snps mode is selected. This will phase bam file without indel calling.', default=False, action='store_true')
    
    phase_group.add_argument("-enable_whatshap",  "--enable_whatshap",  help="Allow WhatsHap to change SNP genotypes when phasing using --distrust-genotypes and --include-homozygous flags (this is not the same as regenotyping), considerably increasing the time needed for phasing. It has a negligible effect on SNP calling accuracy for Nanopore reads, but may make a small improvement for PacBio reads. By default WhatsHap will only phase SNP calls produced by NanoCaller, but not change their genotypes.",  default=False, action='store_true')
    
    args = parser.parse_args()
    
    args.supplementary=False
    
    set_flags=[]
    
    for x in sys.argv:
        if '--' in x:
            set_flags.append(x.replace('-',''))
        
        elif '-' in x:
            set_flags.append(flag_map(x.replace('-','')))
            
    if args.preset:
        for p in preset_dict[args.preset]:
            if p not in set_flags:
                vars(args)[p]=preset_dict[args.preset][p]    
    
    if not args.output:
        args.output=os.getcwd()
    
    os.makedirs(args.output, exist_ok=True)  
    
    
    with open(os.path.join(args.output,'args'),'w') as file:
        file.write('Command: python %s\n\n\n' %(' '.join(sys.argv)))
        file.write('------Parameters Used For Variant Calling------\n')
        for k in vars(args):
            file.write('{}: {}\n'.format(k,vars(args)[k]) )
        
    print('\n%s: Starting NanoCaller.\n\nNanoCaller command and arguments are saved in the following file: %s\n' %(str(datetime.datetime.now()), os.path.join(args.output,'args')), flush=True)
    
    run(args)
    
    elapsed=time.time()-t
    print ('\n%s: Total Time Elapsed: %.2f seconds' %(str(datetime.datetime.now()), elapsed))
