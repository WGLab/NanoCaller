import sys,pysam, time,os,re,copy,argparse,gzip,itertools
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree


def repeat_check(x):
    if len(x)==0:
        return False
    count=1
    letter=x[0]
    i=1
    while count<4 and i<len(x):
        
        if x[i]==letter:
            count+=1
        else:
            letter=x[i]
            count=1
        i+=1
    return count>=4

chrom_length={'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555, 'chr5':181538259, 'chr6':170805979, \
             'chr7':159345973, 'chr8':145138636, 'chr9':138394717, 'chr10':133797422, 'chr11':135086622, 'chr12':133275309,\
             'chr13':114364328, 'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,\
             'chr19':58617616, 'chr20':64444167, 'chr21':46709983, 'chr22':50818468, 'chrX':156040895, 'chrY':57227415}

vcf_path=str(sys.argv[2])
fasta_path='/home/ahsanm/dt_Nanovar/GRCh38.fa'

fastafile=pysam.FastaFile(fasta_path)

chrom=str(sys.argv[1])
start=1
end=chrom_length[chrom]//2

ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }

ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
ref_df.set_index('pos',drop=False,inplace=True)

bcf_in = VariantFile(vcf_path)  # auto-detect input format
gt_map={(0,0):0, (1,1):0, (2,2):0, (1,2):1, (2,1):1, (0,1):1, (1,0):1, (0,2):1,(2,0):1}
indel_pos={}
for rec in bcf_in.fetch(chrom,start,end+1):
    if rec.samples.items()[0][1].get('GT') in gt_map.keys():
        indel_pos[rec.pos]=(rec.samples.items()[0][1].get('GT'),rec.alleles[0],rec.alleles[1])

        
start=chrom_length[chrom]//2
end=chrom_length[chrom]
ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }

ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
ref_df.set_index('pos',drop=False,inplace=True)

for rec in bcf_in.fetch(chrom,start,end+1):
    if rec.samples.items()[0][1].get('GT') in gt_map.keys():
        indel_pos[rec.pos]=(rec.samples.items()[0][1].get('GT'),rec.alleles[0],rec.alleles[1])
        
        
with open(os.path.join(str(sys.argv[3]),'%s.rem_list'%chrom),'w') as file:
    for v_pos in indel_pos.keys():
        if len(indel_pos[v_pos][1])==1 and len(indel_pos[v_pos][1])==1:
            file.write('%s\t%d\n' %(chrom,v_pos))
        elif len(indel_pos[v_pos][1])< len(indel_pos[v_pos][1]):
            seq=''.join(ref_df.reindex(range(v_pos-2,v_pos),axis=0).ref.to_list())+indel_pos[v_pos][-1]+''.join(ref_df.reindex(range(v_pos+1,v_pos+4),axis=0).ref.to_list())
            if repeat_check(seq):
                file.write('%s\t%d\n' %(chrom,v_pos))
                
        elif len(indel_pos[v_pos][1])> len(indel_pos[v_pos][1]):
            seq=''.join(ref_df.reindex(range(v_pos-2,v_pos+4),axis=0).ref.to_list())
            if repeat_check(seq):
                file.write('%s\t%d\n' %(chrom,v_pos))
