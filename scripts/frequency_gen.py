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
    
    bed_path=dct['bed']
    bed_file=os.path.join(bed_path,'%s.bed' %chrom)
    with open(bed_file) as file:
        content=[x.rstrip('\n') for x in file]

    content=[x.split('\t')[1:] for x in content]
    content=[(int(x[0]),int(x[1])) for x in content]
    t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    output=[]
    for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            v_pos=pcol.pos+1
            if t[v_pos]:
                n=pcol.get_num_aligned()
                r=ref_dict[pcol.pos+1]
                if r in 'AGTC' and n>=4:
                    seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                    name=pcol.get_query_names()

                    hap0_bases=[x for x,y in zip(seq,name) if y in dct['hap0']]
                    hap1_bases=[x for x,y in zip(seq,name) if y in dct['hap1']]

                    l1,l2,l3=len(seq),len(hap0_bases),len(hap1_bases)
                    hap0_bases=''.join(hap0_bases).replace('*','')
                    hap1_bases=''.join(hap1_bases).replace('*','')
                    seq=seq.replace('*','')

                    tot_freq=max([x[1] for x in Counter(seq).items() if x[0]!=r]+[0])/len(seq) if len(seq)>0 else 0

                    hap0_freq=max([x[1] for x in Counter(hap0_bases).items() if x[0]!=r]+[0])/len(hap0_bases) if len(hap0_bases)>0 else 0
                    hap1_freq=max([x[1] for x in Counter(hap1_bases).items() if x[0]!=r]+[0])/len(hap1_bases) if len(hap1_bases)>0 else 0

                    if hap0_freq>=0.1 or hap1_freq>=0.1:
                        output.append((v_pos,tot_freq,hap0_freq,hap1_freq,l1,l2,l3))
    return output

def get_5_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    
    samfile = pysam.Samfile(dct['sam_path'], "rb")

    bcf_in = VariantFile(dct['vcf_path'])
    tr_pos={}
    for rec in bcf_in.fetch(chrom,start,end+1):
        tr_pos[rec.pos]=rec.alleles[1]
    
    output=[]
    
    if len(tr_pos)==0:
        return output
    
    for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
        
    

        v_pos=pcol.pos+1

        if v_pos in tr_pos.keys():
            seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
            name=pcol.get_query_names()

            hap0_bases=[x for x,y in zip(seq,name) if y in dct['hap0']]
            hap1_bases=[x for x,y in zip(seq,name) if y in dct['hap1']]

            hap0_bases=''.join(hap0_bases).replace('*','')
            hap1_bases=''.join(hap1_bases).replace('*','')
            seq=seq.replace('*','')

            alt=tr_pos[v_pos]
            tot_freq=seq.count(alt)/len(seq) if len(seq)>0 else 0

            hap0_freq=hap0_bases.count(alt)/len(hap0_bases) if len(hap0_bases)>0 else 0
            hap1_freq=hap1_bases.count(alt)/len(hap1_bases) if len(hap1_bases)>0 else 0

            output.append((v_pos,tot_freq,hap0_freq,hap1_freq))

    return output        

def generate(params,mode='training'):
    cores=params['cpu']
    mode=params['mode']
    print('starting pileups',flush=True)
    pool = mp.Pool(processes=cores)
    fname='%s.pileups.test.freq.phased' %(params['chrom'])
    fname=os.path.join(params['out_path'],fname)
    start,end=params['start'],params['end']
    print('generating frequencies for region %s:%d-%d'%(chrom,start,end),flush=True)
    
    df=pd.read_csv('/home/ahsanm/dt_Nanovar/ground_truth/pedigree/HG001/HG001.%s.bam' %params['chrom'],delim_whitespace=True)
    
    hap0=set(df[(df.haplotype==0)]['#readname'])
    hap1=set(df[(df.haplotype==1)]['#readname'])
    
    params['hap0']=hap0
    params['hap1']=hap1
    
    with open(fname , "w") as file:
        file.write('pos,tot_freq,hap0_freq,hap1_freq,len_seq,len_hap0,len_hap1\n')
        for mbase in range(start,end,int(4e7)):
            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(4e7)),100000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(end,k+100000)
                in_dict_list.append(d)
            results = pool.map(get_training_candidates, in_dict_list)

            
            for res_tuple in results:
                if res_tuple:
                    for stuff in res_tuple:#[pcol.pos+1,n,allele[0],allele[1],alt_freq])
                        file.write('%d,%.3f,%.3f,%.3f,%d,%d,%d\n' %(stuff[0],stuff[1],stuff[2],stuff[3],stuff[4],stuff[5],stuff[6]))


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
    

    
