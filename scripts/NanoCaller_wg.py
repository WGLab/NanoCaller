import time, argparse, os, shutil, sys, NanoCaller
import multiprocessing as mp

if __name__ == '__main__':
    t=time.time()
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('Required arguments')
    parser.add_argument("-mode",  "--mode",  help="Testing mode, options are 'snps', 'indels' and 'both'", type=str, default='both')
    parser.add_argument("-seq",  "--sequencing",  help="Sequencing type, options are 'ont' and 'pacbio'", type=str, default='ont')
    parser.add_argument("-model",  "--model",  help="NanoCaller SNP model to be used, options are 'NanoCaller1' (trained on HG001 Nanopore reads), 'NanoCaller2' (trained on HG002 Nanopore reads) and 'NanoCaller3' (trained on HG003 PacBio reads) ", default='NanoCaller1')
    parser.add_argument("-vcf",  "--vcf",  help="VCF output path, default is current working directory", type=str)
    
    parser.add_argument("-cpu",  "--cpu",  help="CPUs", type=int, default=1)
    parser.add_argument("-min_allele_freq",  "--min_allele_freq",  help="minimum alternative allele frequency", type=float,  default=0.15)
    parser.add_argument("-min_nbr_sites",  "--min_nbr_sites",  help="minimum number of nbr sites", type=int,  default =1)
    
    requiredNamed.add_argument("-bam",  "--bam",  help="Bam file, should be phased if 'indel' mode is selected", required=True)
    requiredNamed.add_argument("-ref",  "--ref",  help="reference genome file with .fai index", required=True)
    
    requiredNamed.add_argument("-prefix",  "--prefix",  help="VCF file prefix", type=str, required=True)
    parser.add_argument("-sample",  "--sample",  help="VCF file sample name", type=str, default='SAMPLE')
    parser.add_argument("-sup",  "--supplementary",  help="Use supplementary reads", type=bool, default=False)
    parser.add_argument("-mincov",  "--mincov",  help="min coverage", type=int, default=8)
    parser.add_argument("-maxcov",  "--maxcov",  help="max coverage", type=int, default=160)
    
    parser.add_argument("-start",  "--start",  help="start, default is 1", type=int)
    parser.add_argument("-end",  "--end",  help="end, default is the end of contig", type=int)
    parser.add_argument("-nbr_t",  "--neighbor_threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended and for PacBio reads '0.3,0.7' is recommended", type=str, default='0.4,0.6')
    parser.add_argument("-ins_t", "--ins_threshold", help="Insertion Threshold",type=float,default=0.4)
    parser.add_argument("-del_t", "--del_threshold", help="Deletion Threshold",type=float,default=0.6)
    
    parser.add_argument("-allow_whatshap",  "--allow_whatshap",  help="Allow WhatsHap to change SNP genotypes when phasing", type=bool,  default=True)
    
    args = parser.parse_args()

    run(args)
    
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)