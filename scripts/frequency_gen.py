import sys,pysam, time,os,re,copy,argparse,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from matplotlib import pyplot as plt
from intervaltree import Interval, IntervalTree


def get_training_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']

    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)
    
    bed_path='/home/ahsanm/dt_Nanovar/bed_by_chrom/HG002/hg38/'
    bed_file=bed_path+'chr%d.bed' %chrom
    with open(bed_file) as file:
        content=[x.rstrip('\n') for x in file]

    content=[x.split('\t')[1:] for x in content]
    content=[(int(x[0]),int(x[1])) for x in content]
    t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)

    rlist=[s for s in fastafile.fetch(chrom,start-1,end-1)]
    output=[]
    for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            n=pcol.get_num_aligned()
            r=rlist[pcol.pos+1-start]

            if r!='N' and n>=4:
                try:
                    seq=''.join(pcol.get_query_sequences()).upper()

                except AssertionError:
                    seq=''.join([pread.alignment.query_sequence[pread.query_position] if pread.query_position else '' for pread in pcol.pileups ])
                alt_freq=max([x[1] for x in Counter(seq).items() if x[0]!=r]+[0])/n
                
                if alt_freq>0.05:
                    output.append([pcol.pos,n,alt_freq])
    return output

def get_5_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']

    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    bed_path='/home/ahsanm/dt_Nanovar/bed_by_chrom/HG002/hg38/'
    bed_file=bed_path+'%s.bed' %chrom
    with open(bed_file) as file:
        content=[x.rstrip('\n') for x in file]

    content=[x.split('\t')[1:] for x in content]
    content=[(int(x[0]),int(x[1])) for x in content]
    t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)

    rlist=[s for s in fastafile.fetch(chrom,start-1,end-1)]
    output=[]
    for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            if bed_path:
                    if not t[pcol.pos+1]:
                        continue
                        
            n=pcol.get_num_aligned()
            r=rlist[pcol.pos+1-start]
            
            if r!='N':
                seq=''.join(pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=False)).upper()
                c=Counter(seq)
                cc=sorted(c, key=c.get)
                try:
                    cc.remove(r)
                except ValueError:
                    pass
                if cc:
                    allele=(cc[-1], c[cc[-1]])

                    alt_freq=allele[1]/n

                    if allele[1]>1 and alt_freq>0.05 and n>=4:
                        output.append([pcol.pos+1,n,allele[0],allele[1],alt_freq])
    return output        

def generate(params,mode='training'):
    cores=params['cpu']
    mode=params['mode']
    print('starting pileups',flush=True)
    pool = mp.Pool(processes=cores)
    fname='%s.pileups.freq_v3' %(params['chrom'])
    fname=os.path.join(params['out_path'],fname)
    start,end=params['start'],params['end']
    print('generating frequencies for region chr%s:%d-%d'%(chrom,start,end),flush=True)
    with open(fname , "w") as file:
        for mbase in range(start,end,int(4e7)):
            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(4e7)),100000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(end,k+100000)
                in_dict_list.append(d)
            results = pool.map(get_5_candidates, in_dict_list)

            for res_tuple in results:
                if res_tuple:
                    for stuff in res_tuple:#[pcol.pos+1,n,allele[0],allele[1],alt_freq])
                        file.write('%d,%d,%s,%d,%.3f\n' %(stuff[0],stuff[1],stuff[2],stuff[3],stuff[4]))


if __name__ == '__main__':
    chrom_length={'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555, 'chr5':181538259, 'chr6':170805979, \
             'chr7':159345973, 'chr8':145138636, 'chr9':138394717, 'chr10':133797422, 'chr11':135086622, 'chr12':133275309,\
             'chr13':114364328, 'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,\
             'chr19':58617616, 'chr20':64444167, 'chr21':46709983, 'chr22':50818468, 'chrX':156040895, 'chrY':57227415}
    parser = argparse.ArgumentParser()

    #-r chromosome region   -m mode   -bam bam file   -ref reference file   -vcf ground truth variants   -o output path
    parser.add_argument("-r", "--region", help="Chromosome region")
    parser.add_argument("-m", "--mode", help="Mode")
    parser.add_argument("-bam", "--bam", help="Bam file")
    parser.add_argument("-ref", "--ref", help="Size")
    parser.add_argument("-vcf", "--vcf", help="Ground truth variants")
    parser.add_argument("-o", "--output", help="Output path")
    parser.add_argument("-w", "--window", help="Window",type=int)
    parser.add_argument("-cpu", "--cpu", help="CPUs",type=int)
    parser.add_argument("-d", "--depth", help="Depth",type=int)
    parser.add_argument("-t", "--threshold", help="Threshold",type=float)
    parser.add_argument("-bed", "--bed", help="BED file")
    parser.add_argument("-mincov", "--mincov", help="min coverage",type=int)
    
    
    args = parser.parse_args()
    
    if len(args.region.split(':'))==2:
        chrom,region=args.region.split(':')
        start,end=int(region.split('-')[0]),int(region.split('-')[1])
        
    else:
        chrom=args.region.split(':')[0]
        start,end=1,chrom_length[chrom]
        
    in_dict={'mode':args.mode, 'chrom':chrom,'start':start,'end':end,\
         'sam_path':args.bam, 'fasta_path':args.ref, 'vcf_path':args.vcf,\
             'out_path':args.output, 'window':args.window, 'depth':args.depth,\
             'threshold':args.threshold, 'cpu':args.cpu, 'bed':args.bed,\
            'mincov':args.mincov}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
    

    
