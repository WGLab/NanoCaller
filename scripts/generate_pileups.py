import sys,pysam, time,os,re,copy,argparse,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from matplotlib import pyplot as plt

def extract_vcf(dct):
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    vcf_path=dct['vcf_path']
    fasta_path=dct['fasta_path']

    
    mapping={'A':0,'G':1,'T':2,'C':3}
    bcf_in = VariantFile(vcf_path)  # auto-detect input format
    fastafile=pysam.FastaFile(fasta_path)
    for rec in bcf_in.fetch():
        prfx='chr' if 'chr' in rec.contig else ''
        chrom=re.findall("(.*?)\d",rec.contig)[0]+re.findall("\d+",chrom)[0]
        break


    d={}
    for rec in bcf_in.fetch(chrom,start,end):
        if len(''.join(rec.alleles))==len(rec.alleles):
            alleles=rec.samples.items()[0][1].items()[0][1]
            d[rec.pos]=[alleles[0]!=alleles[1], mapping[rec.alleles[1]],mapping[rec.ref]]

    df=pd.DataFrame.from_dict(d,orient='index')
    df.reset_index(inplace=True)
    df.rename(columns={'index':'POS',0:'gtype',1:'ALL',2:'REF'},inplace=True)
    #df.sort_values('POS',inplace=True)
    return df

def get_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    threshold=dct['threshold']

    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    fasta_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    sam_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    rlist=[s for s in fastafile.fetch(fasta_chr,start-1,end-1)]
    output=[]
    for pcol in samfile.pileup(sam_chr,start-1,end-1,min_base_quality=0,stepper='nofilter',\
                                           flag_filter=0x800,truncate=True):
            alt_freq=None
            
            try:
                n=pcol.get_num_aligned()
                r=rlist[pcol.pos+1-start]
                
                if r!='N' and n>=16:
                    r=rlist[pcol.pos+1-start]
                    seq=''.join([x for x in pcol.get_query_sequences() if x not in ['','N','n']]).upper()
                    alt_freq=[x[1] for x in Counter(seq).items() if (x[1]>=n*threshold and x[0]!=r)]
                    if alt_freq:
                        output.append((pcol.pos+1,r,r))
                        
            except AssertionError:
                continue
    if output:
        candidates_df=pd.DataFrame(output)
        candidates_df.rename(columns={0:'POS',1:'ALL',2:'REF'},inplace=True)
        mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
        candidates_df=candidates_df[candidates_df['ALL'].map(lambda x: x in ['A','G','T','C'])] 
        
        candidates_df[['ALL','REF']]=candidates_df[['ALL','REF']].applymap(lambda x:mapping[x])
        #candidates_df.rename(columns={0:'POS',1:'ALT_FREQ',2:'REF_FREQ',3:'DEPTH',4:'ALT_NUM',5:'REF_NUM'},inplace=True)
        return candidates_df
    else:
        return pd.DataFrame()
    
def create_training_pileup(in_dct):
    tr_df=extract_vcf(in_dct)
    c_df=get_candidates(in_dct)
    
    tr_list=[]
    tr_pos=[]
    if tr_df.shape[0]!=0:

        tr_df['gtype']=tr_df['gtype'].astype(int)
        tr_pos=list(tr_df.POS)
        tr_df.set_index('POS',drop=False,inplace=True)
        
        for v_pos in tr_df.POS:
            res=create_pileup(in_dct,v_pos=v_pos)
            if res:
                tr_list.append((v_pos,tr_df['gtype'].loc[v_pos],tr_df['ALL'].loc[v_pos],\
                          tr_df['REF'].loc[v_pos],res[0]))


    

    neg_list=[]
    
    if c_df.shape[0]!=0:
        c_df.set_index('POS',drop=False,inplace=True)
        c_df.drop([x for x in tr_pos if x in c_df.index],inplace=True)
        np.random.seed(42)
        sample=np.random.choice(c_df.shape[0], min(c_df.shape[0],5*len(tr_pos)), replace=False)
        c_df=c_df.iloc[sample]
        c_df['gtype']=0
        
        for v_pos in c_df.POS:

            res=create_pileup(in_dct,v_pos=v_pos)
            if res:
                neg_list.append((v_pos,c_df['gtype'].loc[v_pos],c_df['ALL'].loc[v_pos],\
                      c_df['REF'].loc[v_pos],res[0]))
    #gtype, matrix, reads

    return (tr_list,neg_list)

def create_testing_pileup(in_dct):
    test_sites_df=get_candidates(in_dct)
    if test_sites_df.shape[0]==0:
        return None

    test_sites_df.set_index('POS',drop=False,inplace=True)
    pileup_list=[]

    d={}
    for v_pos in test_sites_df.POS:
        
        res=create_pileup(in_dct,v_pos=v_pos)
        if res:
            pileup_list.append(res[0])
            d[v_pos]=(test_sites_df['REF'].loc[v_pos],res[0])
    
    return d
        
def create_pileup(dct,v_pos=None):
    chrom=dct['chrom']
    
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)
    
    #this function creates pileup dataframe
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
    fasta_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    sam_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    
    try:
        du_lim=dct['depth']
        window=dct['window']
        
        v_start=v_pos-window
        v_end=v_pos+window
        
        d={x:{} for x in range(v_start,v_end)}
        features={x:{} for x in range(v_start,v_end)}
        try:
            rlist=[mapping[s] for s in fastafile.fetch(fasta_chr,v_start-1,v_end)]
            
        except KeyError:
            return None
        for pcol in samfile.pileup(chrom,v_start-1,v_end,min_base_quality=0,stepper='nofilter',\
                                           flag_filter=0x800,truncate=True):
            name=pcol.get_query_names()
            seq=pcol.get_query_sequences()
            qual=pcol.get_query_qualities()

            d[pcol.pos+1]={n:s.upper() if len(s)>0 else '*' for (n,s) in zip(name,seq)}
            features[pcol.pos+1]={n:q for (n,q) in zip(name,qual)}


        p_df=pd.DataFrame.from_dict(d)
        p_df.dropna(subset=[v_pos],inplace=True)
        
        if p_df.shape[0]>du_lim:
            p_df=p_df.sample(n=du_lim,replace=False, random_state=1)
        p_df.sort_values(v_pos,inplace=True,ascending=False)
    
        p_df.fillna('N',inplace=True)
        p_df=p_df.applymap(lambda x: mapping[x])
        p_mat=np.array(p_df)
        rnames=list(p_df.index)

        #generatre quality df
        f_df=pd.DataFrame.from_dict(features)
        f_df.fillna(0,inplace=True)
        f_df=f_df.reindex(p_df.index)

        ref_match=(p_mat==rlist)

        #tmp2=2*np.array(ref_match)[:,:,np.newaxis]-1
        tmp=np.dstack([(p_mat==i) for i in range(5)])
        
        #data=np.multiply(tmp,tmp2).astype(np.int8)
        try:
            data=np.dstack((tmp,ref_match,np.array(f_df))).astype(np.int8)
        except ValueError:
            return None
        
        if data.shape[0]<du_lim:
            tmp=np.zeros((du_lim-data.shape[0],data.shape[1],data.shape[2]))
            data=np.vstack((data,tmp))
            rnames+=['Buff' for i in range(tmp.shape[0])]
            
        return (data,rnames)
    except AssertionError:
        return None
        
def generate(params,mode='training'):
    mode=params['mode']
    if mode=='training':
        print('starting pileups')
        cores=mp.cpu_count()
        pool = mp.Pool(processes=mp.cpu_count())
        fname='%s_%d_%d_pileups' %(params['chrom'],params['start'],params['end'])
        true_fname=os.path.join(params['out_path'],fname+'_pos.gz')
        neg_fname=os.path.join(params['out_path'],fname+'_neg.gz')
        p_file=gzip.open(true_fname , "wb")
        n_file=gzip.open(neg_fname , "wb")
        start,end=params['start'],params['end']
        
        
        
        for mbase in range(start,end,int(1e6)):
            print('starting pool:'+str(mbase))
            t=time.time()                

            in_dict_list=[]
            chunk_size=min(100000,int(1e6/cores))
            for k in range(mbase,min(end,mbase+int(1e6)),chunk_size):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=k+chunk_size
                in_dict_list.append(d)
            results_dict = pool.map(create_training_pileup, in_dict_list)
            for result in results_dict:
                if result[1] and (params['examples']=='neg' or params['examples']=='all'):
                    for data in result[1]:
                            pos,gt,allele,ref,mat=data

                            ft=mat[:,:,-1].reshape(-1)
                            mat=mat[:,:,:-1].reshape(-1)
                            s='%d:%d:%d:%d:%s:%s\n' %(pos,gt,allele,ref,''.join(mat.astype('<U1')),'|'.join(ft.astype('<U2')))
                            n_file.write(s.encode('utf-8'))

                if result[0] and (params['examples']=='pos' or params['examples']=='all'):
                    for data in result[0]:
                        pos,gt,allele,ref,mat=data

                        ft=mat[:,:,-1].reshape(-1)
                        mat=mat[:,:,:-1].reshape(-1)
                        s='%d:%d:%d:%d:%s:%s\n' %(pos,gt,allele,ref,''.join(mat.astype('<U1')),'|'.join(ft.astype('<U2')))
                        p_file.write(s.encode('utf-8'))

            
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed)
            print('finishing pool:'+str(mbase))
            
    elif mode=='testing':
        print('starting pileups')
        pool = mp.Pool(processes=mp.cpu_count())
        fname='%s_%d_%d_pileups' %(params['chrom'],params['start'],params['end'])
        fname=os.path.join(params['out_path'],fname+'_multi_test_v2.gz')
        file=gzip.open(fname , "wb")
        start,end=params['start'],params['end']
        for mbase in range(start,end,int(1e4)):
            print('starting pool:'+str(mbase))
            t=time.time()                

            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(1e4)),1000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=k+1000
                in_dict_list.append(d)
            results = pool.map(create_testing_pileup, in_dict_list)
            for res_dict in results:
                if res_dict:
                    for k,p in res_dict.items():
                        ref,mat=p

                        ft=mat[:,:,-1].reshape(-1)
                        mat=mat[:,:,:-1].reshape(-1)
                        s='%d:%d:%s:%s\n' %(k,ref, ''.join(mat.astype('<U1')),'|'.join(ft.astype('<U2')))
                        file.write(s.encode('utf-8'))

            
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed)
            print('finishing pool:'+str(mbase))



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    #-r chromosome region   -m mode   -bam bam file   -ref reference file   -vcf ground truth variants   -o output path
    parser.add_argument("-r", "--region", help="Chromosome region")
    parser.add_argument("-m", "--mode", help="Mode")
    parser.add_argument("-bam", "--bam", help="Bam file")
    parser.add_argument("-ref", "--ref", help="Size")
    parser.add_argument("-vcf", "--vcf", help="Ground truth variants")
    parser.add_argument("-o", "--output", help="Output path")
    parser.add_argument("-w", "--window", help="Window",type=int)
    parser.add_argument("-d", "--depth", help="Depth",type=int)
    parser.add_argument("-t", "--threshold", help="Threshold",type=float)
    parser.add_argument("-ex", "--examples", help="Threshold",type=str)
    
    args = parser.parse_args()
    
    chrom,region=args.region.split(':')
    start,end=int(region.split('-')[0]),int(region.split('-')[1])
    
    in_dict={'mode':args.mode, 'chrom':chrom,'start':start,'end':end,\
         'sam_path':args.bam,'fasta_path':args.ref,'vcf_path':args.vcf,'out_path':args.output,'window':args.window,'depth':args.depth,'threshold':args.threshold,'examples':args.examples}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
    
