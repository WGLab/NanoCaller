import sys,pysam, time,os,re,copy,argparse,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from matplotlib import pyplot as plt
from intervaltree import Interval, IntervalTree


def get_nbr(dct):
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    vcf_path=dct['vcf_path']
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    
    if dct['type'] in ['training','train']:
        bcf_in = VariantFile(vcf_path)  # auto-detect input format
        fastafile=pysam.FastaFile(fasta_path)

        gt_map={(0,0):0, (1,1):0, (2,2):0, (1,2):1, (2,1):1, (0,1):1, (1,0):1, (0,2):1,(2,0):1,(1,None):0}
        tr_pos={}
        for rec in bcf_in.fetch(chrom,start,end+1):
            if gt_map[rec.samples.items()[0][1].get('GT')]:
                tr_pos[rec.pos]=rec.ref
            
        output=[]
        for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
                if pcol.pos+1 in tr_pos.keys():
                    n=pcol.get_num_aligned()
                    r=tr_pos[pcol.pos+1]

                    if r!='N' and n>=dct['mincov']:
                        seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                        name=pcol.get_query_names()
                        output.append([pcol.pos+1,r,seq,':'.join(name)])
        return output
    
    elif dct['type'] in ['test','testing']:
        rlist=[s for s in fastafile.fetch(chrom,start-1,end-1)]
        output=[]
        for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
                n=pcol.get_num_aligned()
                r=rlist[pcol.pos+1-start]

                if r!='N' and n>=dct['mincov']:
                    seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                    alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=r and x[0] in 'AGTC')]+[0])/n

                    if 0.3<=alt_freq and alt_freq<0.7:
                        name=pcol.get_query_names()
                        output.append([pcol.pos+1,n,alt_freq,r,seq,':'.join(name)])
        return output
    
    elif dct['type']=='redo':
        output=[]
        bcf_in = VariantFile(dct['vcf_path'])
        cand={}
        for rec in bcf_in.fetch(chrom,start,end):
            if rec.samples.items()[0][1].get('GT')==(0,1) and rec.samples.items()[0][1].get('GP')>=5:
                cand[rec.pos]=(rec.ref,'0/1')
            elif rec.samples.items()[0][1].get('GT') == (1,0):
                cand[rec.pos]=(rec.ref,'1/0')
                
        for pcol in samfile.pileup(chrom,start-1,end,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            if pcol.pos+1 in cand.keys():
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()
                output.append([pcol.pos+1,cand[pcol.pos+1][0],cand[pcol.pos+1][1],seq,':'.join(name)])
        return output

            
        
def generate(params,mode='training'):
    cores=params['cpu']
    print('starting pileups',flush=True)
    pool = mp.Pool(processes=cores)
    
    start,end=params['start'],params['end']
    if start==None:
        start=1
    if end==None:
        end=chrom_length[params['chrom']]
        
    print('generating frequencies for %s'%(chrom),flush=True)
    
    if params['type'] in ['training','train']:
        fname='%s.pileups.neighbors.train' %params['chrom']
        fname=os.path.join(params['out_path'],fname)
        with open(fname , "w") as file:
            file.write('pos,ref,seq,names\n')
            for mbase in range(start,end,int(4e6)):
                in_dict_list=[]
                for k in range(mbase,min(end,mbase+int(4e6)),100000):
                    d = copy.deepcopy(params)
                    d['start']=k
                    d['end']=min(end,k+100000)
                    in_dict_list.append(d)
                results = pool.map(get_nbr, in_dict_list)

                for res_tuple in results:
                    if res_tuple:
                        for stuff in res_tuple:
                            file.write('%d,%s,%s,%s\n' %(stuff[0],stuff[1],stuff[2],stuff[3]))

                           
    elif params['type']  in ['test','testing']:
        fname='%s.pileups.neighbors.test' %params['chrom']
        fname=os.path.join(params['out_path'],fname)
        with open(fname , "w") as file:
            file.write('pos,depth,freq,ref,seq,names\n')
            for mbase in range(start,end,int(4e6)):
                in_dict_list=[]
                for k in range(mbase,min(end,mbase+int(4e6)),100000):
                    d = copy.deepcopy(params)
                    d['start']=k
                    d['end']=min(end,k+100000)
                    in_dict_list.append(d)
                results = pool.map(get_nbr, in_dict_list)
                for res in results:
                    if res:
                        for stuff in res:
                            file.write('%d,%d,%.3f,%s,%s,%s\n' %(stuff[0],stuff[1],stuff[2],stuff[3],stuff[4],stuff[5]))
    
    elif params['type']=='redo':
        fname='%s.pileups.neighbors.%s' %(params['chrom'],params['suffix'])
        fname=os.path.join(params['out_path'],fname)
        with open(fname , "w") as file:
                    file.write('pos,ref,gt,seq,names\n')
                    for mbase in range(start,end,int(4e6)):
                        in_dict_list=[]
                        for k in range(mbase,min(end,mbase+int(4e6)),100000):
                            d = copy.deepcopy(params)
                            d['start']=k
                            d['end']=min(end,k+100000)
                            d['type']='redo'
                            in_dict_list.append(d)
                        results = pool.map(get_nbr, in_dict_list)
                        for res in results:
                            if res:
                                for stuff in res:
                                    file.write('%d,%s,%s,%s,%s\n' %(stuff[0],stuff[1],stuff[2],stuff[3],stuff[4]))
        return fname
                                
if __name__ == '__main__':
    chrom_length={'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555, 'chr5':181538259, 'chr6':170805979, \
             'chr7':159345973, 'chr8':145138636, 'chr9':138394717, 'chr10':133797422, 'chr11':135086622, 'chr12':133275309,\
             'chr13':114364328, 'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,\
             'chr19':58617616, 'chr20':64444167, 'chr21':46709983, 'chr22':50818468, 'chrX':156040895, 'chrY':57227415}
    parser = argparse.ArgumentParser()

    #-r chromosome region   -m mode   -bam bam file   -ref reference file   -vcf ground truth variants   -o output path
    parser.add_argument("-r", "--region", help="Chromosome region")
    parser.add_argument("-bam", "--bam", help="Bam file")
    parser.add_argument("-ref", "--ref", help="Size")
    parser.add_argument("-o", "--output", help="Output path")
    parser.add_argument("-w", "--window", help="Window",type=int)
    parser.add_argument("-cpu", "--cpu", help="CPUs",type=int)
    parser.add_argument("-d", "--depth", help="Depth",type=int)
    parser.add_argument("-t", "--threshold", help="Threshold",type=float)
    parser.add_argument("-bed", "--bed", help="BED file")
    parser.add_argument("-mincov", "--mincov", help="min coverage",type=int)
    parser.add_argument("-type", "--type", help="Nbr type")
    parser.add_argument("-vcf", "--vcf", help="Ground truth variants")
    
    args = parser.parse_args()
    
    if len(args.region.split(':'))==2:
        chrom,region=args.region.split(':')
        start,end=int(region.split('-')[0]),int(region.split('-')[1])
        
    else:
        chrom=args.region.split(':')[0]
        start,end=1,chrom_length[chrom]
        
    in_dict={'chrom':chrom,'start':start,'end':end,\
         'sam_path':args.bam, 'fasta_path':args.ref,\
             'out_path':args.output, 'window':args.window, 'depth':args.depth,\
             'threshold':args.threshold, 'cpu':args.cpu, 'bed':args.bed,\
            'mincov':args.mincov,'type':args.type,'vcf_path':args.vcf}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
    

    
