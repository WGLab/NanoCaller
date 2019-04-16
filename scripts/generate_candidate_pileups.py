import pysam, sys
from collections import Counter
import pandas as pd
import numpy as np


def get_candidate(chr_num,start,end,du_lim,dl_lim,threshold,samfile):
    output=[]
    for pcol in samfile.pileup(chr_num,start-1,end,min_base_quality=0,stepper='nofilter',\
                                       flag_filter=0x800,truncate=True):
        
        #if pileupcolumn.n<dl_lim:                
        n=pcol.get_num_aligned()
        xx=[x.upper() for x in pcol.get_query_sequences()]
        if n>=dl_lim and max(Counter(xx).values())/n<=threshold:
            output.append(pcol.pos+1)
    return output

def create_pileup(chr_num,var_pos,du_lim,window,samfile,fastafile):
    
    #this function creates pileup dataframe
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}

    d={x:{} for x in range(var_pos-window,var_pos+window)}
    features={x:{} for x in range(var_pos-window,var_pos+window)}
    rlist=[mapping[s] for s in fastafile.fetch('chr20',var_pos-1-window,var_pos+window-1)]
    for pcol in samfile.pileup(chr_num,var_pos-1-window,var_pos+window-1,min_base_quality=0,stepper='nofilter',\
                                       flag_filter=0x800,truncate=True):
        name=pcol.get_query_names()
        seq=pcol.get_query_sequences()
        qual=pcol.get_query_qualities()

        d[pcol.pos+1]={n:s.upper() if len(s)>0 else '*' for (n,s) in zip(name,seq)}
        features[pcol.pos+1]={n:q for (n,q) in zip(name,qual)}

    #join reference and reads, convert letters to number and then to rgb color channels
    
    p_df=pd.DataFrame.from_dict(d)
    p_df.dropna(subset=[var_pos],inplace=True)
    
    if p_df.shape[0]>du_lim:
        p_df=p_df.sample(n=du_lim,replace=False, random_state=1)
    p_df.sort_values(var_pos,inplace=True,ascending=False)
    
    p_df.fillna('N',inplace=True)
    p_df=p_df.applymap(lambda x: mapping[x])
    p_mat=np.array(p_df)
    
    
    #generatre quality df
    f_df=pd.DataFrame.from_dict(features)
    f_df.dropna(subset=[var_pos],inplace=True)
    f_df.fillna(0,inplace=True)
    f_df=f_df.reindex(p_df.index)

    ref_match=(p_mat==rlist).astype(int)
    
    tmp=np.dstack([(p_mat==i).astype(int) for i in range(5)])
    tmp2=np.zeros(p_mat.shape)
    tmp2[:,window]=1
    
    data=np.dstack((tmp,tmp2,np.array(ref_match),np.array(f_df)))

    if data.shape[0]<du_lim:
        tmp=np.zeros((du_lim-data.shape[0],data.shape[1],data.shape[2]))
        data=np.vstack((data,tmp))
    
    return data.reshape(-1)
    
def generate(chr_num,start,end,sam_path,fasta_path):
    window=16
    du_lim=32
    dl_lim=12
    threshold=0.7
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)
    
    candidates=get_candidate(chr_num,start,end,du_lim,dl_lim,threshold,samfile)
    total=[]

    for var_pos in candidates:
        data=create_pileup(chr_num,var_pos,du_lim,window,samfile,fastafile)
        total.append(data)

	
    total=np.array(total)
    return total, candidates


if __name__ == '__main__':
    if len(sys.argv)!=6:
        print('Wrong number of inputs')
        sys.exit(0)
    
    chr_num,start,end,sam_path,fasta_path=sys.argv[1:6]
    pileups,candidates=generate(chr_num,int(start),int(end),sam_path,fasta_path)
    name=chr_num +'_'+ str(start) +'_'+ str(end)+'.npz'
    np.savez(name,pileups=pileups, candidates=candidates)
    
    
