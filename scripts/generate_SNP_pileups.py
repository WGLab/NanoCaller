import sys, pysam, time, os, re, copy, argparse, gzip, itertools, random
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree

base_to_num_map={'*':4,'A':0,'G':1,'T':2,'C':3,'N':4}

def get_nbr(dct):
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']

    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    rlist=[s for s in fastafile.fetch(chrom,start-1,end-1)]
    output=[]
    for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            n=pcol.get_num_aligned()
            r=rlist[pcol.pos+1-start]

            if r!='N' and n>=dct['mincov']:
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=r and x[0] in 'AGTC')]+[0])/n

                if dct['threshold'][0]<=alt_freq and alt_freq<dct['threshold'][1]:
                    name=pcol.get_query_names()
                    output.append([pcol.pos+1,n,alt_freq,r,seq,':'.join(name)])
    return output


def get_testing_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
        
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    threshold=dct['threshold']
    
    nbr_size=20
    
    cnd_df=dct['cnd_df']
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
    
    pileup_dict={}
    
    output={}
    
    
    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            if ref_df.loc[pcol.pos+1].ref in 'AGTC':
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()
                n=pcol.get_num_aligned()

                if n>=dct['mincov'] and pcol.pos+1>=start and pcol.pos+1<=end:
                    alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=ref_df.loc[pcol.pos+1].ref and x[0] in 'AGTC')]+[0])/n

                    if dct['min_allele_freq']<=alt_freq:
                        pileup_dict[pcol.pos+1]={n:s for (n,s) in zip(name,seq)}
                        output[pcol.pos+1]=(n,alt_freq)
                    
    pileup_list=[]
    
    ref_df[['ref']]=ref_df[['ref']].applymap(lambda x:base_to_num_map[x])
    
    pos_list=output.keys()

    output_pos,output_ref,output_mat,output_dp,output_freq=[],[],[],[],[]
    
    if pos_list:

        for v_pos in pos_list:
            
            if dct['seq']=='ont':
                ls=cnd_df[(abs(cnd_df['pos']-v_pos)<50000)].pos
                ls=sorted(ls)


                ls1_0= [p for p in ls if (p>=v_pos-2000) &  (p<v_pos)][:2]
                ls1_1= [p for p in ls if (p>=v_pos-5000) &  (p<v_pos-2000)][-3:]
                ls1_2= [p for p in ls if (p>=v_pos-10000) & (p<v_pos-5000)][-4:]
                ls1_3= [p for p in ls if (p>=v_pos-20000) & (p<v_pos-10000)][-5:]
                ls1_4= [p for p in ls if                  (p<v_pos-20000)][-6:]


                ls2_0= [p for p in ls if (p>v_pos) & (p<=v_pos+2000)][-2:]
                ls2_1= [p for p in ls if (p>v_pos+2000) & (p<=v_pos+5000)][:3]
                ls2_2= [p for p in ls if (p>v_pos+5000) & (p<=v_pos+10000)][:4]
                ls2_3= [p for p in ls if (p>v_pos+10000) & (p<=v_pos+20000)][:5]
                ls2_4= [p for p in ls if (p>v_pos+20000)][:6]
            
            
                ls_total=sorted(ls1_0+ls1_1+ls1_2+ls1_3+ls1_4 + ls2_0+ls2_1+ls2_2+ls2_3+ls2_4)

            else:
                ls=cnd_df[(abs(cnd_df['pos']-v_pos)<20000)].pos
                ls=sorted(ls)
                
                ls1_0= [p for p in ls if (p>=v_pos-20000) &  (p<v_pos)][-20:]
                ls2_0= [p for p in ls if (p>v_pos) & (p<=v_pos+2000)][:20]
                
                ls_total=sorted(ls1_0+ls2_0)
                
            nbr_dict={}

            rlist=[(v_pos,ref_df.loc[v_pos].ref)]
            
            sample=list(pileup_dict[v_pos].keys())
            
            if len(sample) > dct['maxcov']:
                sample=random.sample(sample,min(len(sample), dct['maxcov']))

            sample=set(sample)
            nbr_dict[v_pos]={x:pileup_dict[v_pos][x] for x in sample}
            
            for nb_pos in ls_total:
                tmp_list=cnd_df.loc[nb_pos].names.split(':')
                if set(tmp_list) & sample:
                        nbr_dict[nb_pos]={n:s for (n,s) in zip(tmp_list,cnd_df.loc[nb_pos].seq) if n in sample}
                        rlist.append((nb_pos, cnd_df.loc[nb_pos].ref))

            rlist=sorted(rlist, key=lambda x:x[0])
            total_rlist=np.array([x[1] for x in rlist])
            
            if len(total_rlist)<dct['min_nbr_sites']:
                continue
            
            p_df=pd.DataFrame.from_dict(nbr_dict)
            p_df = p_df.reindex(sorted(p_df.columns), axis=1)
            p_df.fillna('N',inplace=True)
            p_df=p_df.applymap(lambda x: base_to_num_map[x])

            v_ind=p_df.columns.get_loc(v_pos)
            p_mat=np.array(p_df)
            mat=np.dstack([np.sum(np.eye(5)[p_mat[p_mat[:,v_ind]==i]],axis=0) for i in range(4)]).transpose(2,0,1)[:,:,:4]

            total_ref=np.eye(5)[total_rlist.astype(int)]
            total_ref[:,4]=0
            total_ref=total_ref[np.newaxis,:]
            mat=np.dstack([mat,np.zeros([4,mat.shape[1]])+np.eye(4)[ref_df.loc[v_pos].ref][:,np.newaxis]])
            data=np.vstack([total_ref,np.multiply(mat,1-2*total_ref)])
            data=np.hstack([np.zeros([5,nbr_size-sum(p_df.columns<v_pos),5]),data,np.zeros([5, nbr_size- sum(p_df.columns>v_pos), 5])]).astype(np.int8)
            
            output_pos.append(v_pos)
            output_ref.append(ref_df.loc[v_pos].ref)
            output_mat.append(data)
            output_dp.append(output[v_pos][0])
            output_freq.append(output[v_pos][1])
            
    output_mat=np.array(output_mat).astype(np.float32)    
    output_pos=np.array(output_pos)
    
    output_ref=np.eye(max(4,np.max(output_ref)+1))[np.array(output_ref)].astype(np.int8)
    output_ref=output_ref[:,:4]

    output_dp=np.array(output_dp)
    output_freq=np.array(output_freq)
    
    return (output_pos,output_ref,output_mat,output_dp,output_freq)
    
def generate(params):
    cores=params['cpu']
    chrom=params['chrom']
    threshold=params['threshold']
    start,end=params['start'],params['end']    

    
    print('starting pileup generation for %s:%d-%d' %(chrom,start,end) ,flush=True)
    pool = mp.Pool(processes=cores)
    
    in_dict_list=[]
    for k in range(max(1,start-50000),end+50000,100000):
        d = copy.deepcopy(params)
        d['start']=k
        d['end']=min(end+50000,k+100000)
        in_dict_list.append(d)
        
    nbr_results = pool.map(get_nbr, in_dict_list)
    
    cnd_df=pd.concat([pd.DataFrame(lst) for lst in nbr_results])
    if len(cnd_df)==0:
        return ([],[],[],[],[])
    cnd_df.rename(columns={0:'pos',1:'depth',2:'freq',3:'ref',4:'seq', 5:'names'},inplace=True)
   
    cnd_df.set_index('pos',inplace=True,drop=False)
    cnd_df[['ref']]=cnd_df[['ref']].applymap(lambda x:base_to_num_map[x])
    params['cnd_df']=cnd_df
                                    
    print('neighbor candidates generated',flush =True)
    in_dict_list=[]
    for k in range(start,end,100000):
        d = copy.deepcopy(params)
        d['start']=k
        d['end']=min(end,k+100000)
        in_dict_list.append(d)
    results = pool.map(get_testing_candidates, in_dict_list)

    pos=np.vstack([res[0][:,np.newaxis] for res in results if len(res[0])>0])
    mat=np.vstack([res[2] for res in results if len(res[0])>0])
    
    ref=np.vstack([res[1] for res in results if len(res[0])>0 ])
    dp=np.vstack([res[3][:,np.newaxis] for res in results if len(res[0])>0])

    freq=np.vstack([res[4][:,np.newaxis] for res in results if len(res[0])>0])
    
    return pos,mat,ref,dp,freq