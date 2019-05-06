import sys,pysam, time,os,re,copy,argparse,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from matplotlib import pyplot as plt


def extract_vcf(dct,mode='all',VT=['SNP']):
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    vcf_path=dct['vcf_path']
    
    def gtype(num):
        if num[0]==num[1]:
            if num[0]==0:
                return 0
            else:
                return 1
        else:
            return 2

    

    bcf_in = VariantFile(vcf_path)  # auto-detect input format
    for rec in bcf_in.fetch():
        prfx='chr' if 'chr' in rec.contig else ''
        chrom=re.findall("(.*?)\d",rec.contig)[0]+re.findall("\d+",chrom)[0]
        break

    if mode=='common':
        d=[rec.pos for rec in bcf_in.fetch(chrom,start,end)\
                if list(filter(lambda x: x[0]=='VT',rec.info.items()))[0][1][0] in VT]
        return pd.DataFrame(d,columns=['POS'])
    else:
        d=[(rec.pos,gtype(rec.samples.items()[0][1].items()[0][1]))\
                for rec in bcf_in.fetch(chrom,start,end)\
                if len(''.join(rec.alleles))==len(rec.alleles)]
        df=pd.DataFrame(d)
        df.rename(columns={0:'POS',1:'gtype'},inplace=True)
        return df
    


def get_candidates(dct,df=None):
    
    if dct['mode'] =='testing':
        chrom=dct['chrom']
        start=dct['start']
        end=dct['end']
        sam_path=dct['sam_path']
        fasta_path=dct['fasta_path']
        threshold=dct['threshold']
        
        samfile = pysam.Samfile(sam_path, "rb")
        fastafile=pysam.FastaFile(fasta_path)
        name=chrom+'_'+str(start)+'_'+str(end)+'_'+'candidates'
        
        fasta_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
        sam_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
        rlist=[s for s in fastafile.fetch(fasta_chr,start-1,end-1)]
        output=[]
        for pcol in samfile.pileup(sam_chr,start-1,end-1,min_base_quality=0,stepper='nofilter',\
                                               flag_filter=0x800,truncate=True):
                n=pcol.get_num_aligned()
                r=rlist[pcol.pos+1-start]
                seq=[s.upper() if len(s)>0 else '*' for s in pcol.get_query_sequences()]
                alt_freq=[x[1] for x in Counter(seq).items() if (x[1]>=4 and x[1]>=n*threshold and x[0]!=r)]
                if alt_freq:
                        output.append((pcol.pos+1,max(alt_freq)/n,seq.count(r)/n,n,max(alt_freq),seq.count(r)))

        candidates_df=pd.DataFrame(output)
        candidates_df.rename(columns={0:'POS',1:'ALT_FREQ',2:'REF_FREQ',3:'DEPTH',4:'ALT_NUM',5:'REF_NUM'},inplace=True)

        return candidates_df
    
    elif dct['mode']=='training':
        tmp=pd.DataFrame([i+100 for i in df.POS])
        tmp.rename(columns={0:'POS'},inplace=True)
        tmp_df=pd.merge(df,tmp, on='POS',how='outer').fillna(0)
        tmp_df.sort_values('POS',inplace=True)
        return tmp_df
        


def create_training_pileup(in_dct):
    tmp=extract_vcf(in_dct)
    if tmp.shape[0]==0:
        return None
    train_sites_df=get_candidates(in_dct,tmp)

    train_sites_df.set_index('POS',drop=False,inplace=True)
    d={}
    for v_pos in train_sites_df.POS:
        
        res=create_pileup(in_dct,v_pos=v_pos)
        d[v_pos]=(train_sites_df['gtype'].loc[v_pos],res[0],res[1])
        #gtype, matrix, reads
    
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
    
    if dct['mode']=='training':
        du_lim=32
        window=50
        
        v_start=v_pos-window
        v_end=v_pos+window
        
        d={x:{} for x in range(v_start,v_end)}
        features={x:{} for x in range(v_start,v_end)}
        rlist=[mapping[s] for s in fastafile.fetch(fasta_chr,v_start-1,v_end)]
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
        tmp2=np.zeros(p_mat.shape)
        tmp2[:,window]=1
        #tmp2=2*np.array(ref_match)[:,:,np.newaxis]-1
        tmp=np.dstack([(p_mat==i) for i in range(5)])
        
        #data=np.multiply(tmp,tmp2).astype(np.int8)
        data=np.dstack((tmp,tmp2,ref_match,np.array(f_df))).astype(np.int8)
        
        
        if data.shape[0]<du_lim:
            tmp=np.zeros((du_lim-data.shape[0],data.shape[1],data.shape[2]))
            data=np.vstack((data,tmp))
            rnames+=['Buff' for i in range(tmp.shape[0])]
            
        return (data,rnames)
    else:
        start=dct['start']
        end=dct['end']
        #testing code
        
        
def generate(in_dct,mode='training'):
    if mode=='training':
        print('starting pileups')
        pool = mp.Pool(processes=mp.cpu_count())
        fname='%s_%d_%d_pileups' %(in_dct['chrom'],in_dct['start'],in_dct['end'])
        fname=os.path.join(in_dict['out_path'],fname+'.gz')
        file=gzip.open(fname , "wb")
        start,end=in_dict['start'],in_dict['end']
        for mbase in range(start,end,int(1e6)):
            print('starting pool:'+str(mbase))
            t=time.time()                

            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(1e6)),100000):
                d = copy.deepcopy(in_dct)
                d['start']=k
                d['end']=k+100000
                in_dict_list.append(d)
            results = pool.map(create_training_pileup, in_dict_list)
            for res_dict in results:
                if res_dict:
                    for k,p in res_dict.items():
                        a,b,n=p

                        f=b[:,:,-1].reshape(-1)
                        b=b[:,:,:-1].reshape(-1)
                        s='%d:%s:%s:%s:%s\n' %(k,a,'|'.join(n),''.join(b.astype('<U1')),\
                                             '|'.join(f.astype('<U2')))
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

    args = parser.parse_args()
    
    chrom,region=args.region.split(':')
    start,end=int(region.split('-')[0]),int(region.split('-')[1])
    
    in_dict={'mode':args.mode, 'chrom':chrom,'start':start,'end':end,\
         'sam_path':args.bam,'fasta_path':args.ref,'vcf_path':args.vcf,'out_path':args.output}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
    
