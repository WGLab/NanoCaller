import sys,pysam, time,os,re,copy,argparse,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from matplotlib import pyplot as plt
from intervaltree import Interval, IntervalTree

def extract_vcf(dct):
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    vcf_path=dct['vcf_path']
    fasta_path=dct['fasta_path']

    
    mapping={'A':0,'G':1,'T':2,'C':3}
    bcf_in = VariantFile(vcf_path)  # auto-detect input format
    fastafile=pysam.FastaFile(fasta_path)
    '''for rec in bcf_in.fetch():
        prfx='chr' if 'chr' in rec.contig else ''
        chrom=re.findall("(.*?)\d",rec.contig)[0]+re.findall("\d+",chrom)[0]
        break'''


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



def get_training_candidates(dct,tr_pos):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    threshold=dct['threshold']

    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)
    
    bed_file='/home/ahsanm1/umair_wlab/data/NanoVar_data/bed_by_chrom/%s.bed' %chrom
    with open(bed_file) as file:
        content=[x.rstrip('\n') for x in file]
    
    content=[x.split('\t')[1:] for x in content]
    content=[(int(x[0]),int(x[1])) for x in content]
    t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)
    
    
    
    fasta_chr,sam_chr=chrom,chrom
    #fasta_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    #sam_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    rlist=[s for s in fastafile.fetch(fasta_chr,start-1,end-1)]
    output={0:[],5:[],10:[],15:[],20:[],25:[]}
    for pcol in samfile.pileup(sam_chr,start-1,end-1,min_base_quality=0,stepper='nofilter',\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            alt_freq=None
            if not t[pcol.pos+1]:
                continue
            try:
                n=pcol.get_num_aligned()
                r=rlist[pcol.pos+1-start]
                
                if r!='N' and n>=4:
                    r=rlist[pcol.pos+1-start]
                    seq=''.join([x for x in pcol.get_query_sequences() if x not in ['','N','n']]).upper()
                    alt_freq=max([x[1] for x in Counter(seq).items() if x[0]!=r]+[0])/n
                    #alt_freq=[x[1] for x in Counter(seq).items() if (x[1]>=n*threshold and x[0]!=r)]
                    
                    if alt_freq<0.05:
                        output[0].append((pcol.pos+1,r,r))
                        
                    elif 0.05<=alt_freq<0.10:
                        output[5].append((pcol.pos+1,r,r))
                    elif 0.10<=alt_freq<0.15:
                        output[10].append((pcol.pos+1,r,r))
                    elif 0.15<=alt_freq<0.20:
                        output[15].append((pcol.pos+1,r,r))
                    elif 0.20<=alt_freq<0.25:
                        output[20].append((pcol.pos+1,r,r))
                    elif 0.25<=alt_freq:
                        output[25].append((pcol.pos+1,r,r))
                        
            except AssertionError:
                continue
    np.random.seed(42)
    result={}
    sizes={0:2*len(tr_pos), 5:len(tr_pos)//2,10:len(tr_pos)//2,15:len(tr_pos)//2, 20:len(tr_pos)//2, 25:2*len(tr_pos)}
    for i in [0,5,10,15,20,25]:
        pos_list=output[i]
        candidates_df=pd.DataFrame()
        
        if pos_list:
            candidates_df=pd.DataFrame(pos_list)

            candidates_df.rename(columns={0:'POS',1:'ALL',2:'REF'},inplace=True)

            candidates_df.set_index('POS',drop=False,inplace=True)
            candidates_df.drop([x for x in tr_pos if x in candidates_df.index],inplace=True)

            sample=np.random.choice(len(candidates_df), min(len(candidates_df),sizes[i]), replace=False)
            candidates_df=candidates_df.iloc[sample]

            mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
            candidates_df=candidates_df[candidates_df['ALL'].map(lambda x: x in ['A','G','T','C'])] 
            
            try:
                candidates_df[['ALL','REF']]=candidates_df[['ALL','REF']].applymap(lambda x:mapping[x])
            #candidates_df.rename(columns={0:'POS',1:'ALT_FREQ',2:'REF_FREQ',3:'DEPTH', 4:'ALT_NUM', 5:'REF_NUM'}, inplace=True)
            except KeyError:
                pass
                
        result[i]=candidates_df

    return result

    
def create_training_pileup(in_dct):
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
    neg_list={0:[],5:[],10:[],15:[],20:[],25:[]}
    tr_list=[]
    if in_dct['mode']=='training':
        tr_df=extract_vcf(in_dct)
                
        tr_pos=[]
        if tr_df.shape[0]!=0:

            tr_df['gtype']=tr_df['gtype'].astype(int)
            tr_pos=list(tr_df.POS)
            tr_df.set_index('POS',drop=False,inplace=True)

            for v_pos in tr_df.POS:
                res=create_pileup_image(in_dct,v_pos=v_pos)
                if res:
                    tr_list.append((v_pos,tr_df['gtype'].loc[v_pos],tr_df['ALL'].loc[v_pos],\
                              tr_df['REF'].loc[v_pos],res[0]))



            neg_df_list=get_training_candidates(in_dct,tr_pos)

            

            for i in [0,5,10,15,20,25]:
                c_df=neg_df_list[i]

                if c_df.shape[0]!=0:

                    c_df['gtype']=0

                    for v_pos in c_df.POS:

                        res=create_pileup_image(in_dct,v_pos=v_pos)
                        if res:
                            neg_list[i].append((v_pos,c_df['gtype'].loc[v_pos],c_df['ALL'].loc[v_pos],\
                                  c_df['REF'].loc[v_pos],res[0]))
                #gtype, matrix, reads

        return (tr_list,neg_list)    
    




def get_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    threshold=dct['threshold']

    bed_file='/home/ahsanm1/umair_wlab/data/NanoVar_data/bed_by_chrom/%s.bed' %chrom
    with open(bed_file) as file:
        content=[x.rstrip('\n') for x in file]
    
    content=[x.split('\t')[1:] for x in content]
    content=[(int(x[0]),int(x[1])) for x in content]
    t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)
    
    
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    fasta_chr,sam_chr=chrom,chrom
    #fasta_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    #sam_chr=re.findall("(.*?)\d",fastafile.references[0])[0]+re.findall("\d+",chrom)[0]
    rlist=[s for s in fastafile.fetch(fasta_chr,start-1,end-1)]
    output=[]
    for pcol in samfile.pileup(sam_chr,start-1,end-1,min_base_quality=0,stepper='nofilter',\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            alt_freq=None
            if not t[pcol.pos+1]:
                continue
            try:
                n=pcol.get_num_aligned()
                r=rlist[pcol.pos+1-start]
                
                if r!='N' and n>=4:
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
        
        try:
            candidates_df[['ALL','REF']]=candidates_df[['ALL','REF']].applymap(lambda x:mapping[x])
        except KeyError:
            return pd.DataFrame()
        #candidates_df.rename(columns={0:'POS',1:'ALT_FREQ',2:'REF_FREQ',3:'DEPTH',4:'ALT_NUM',5:'REF_NUM'},inplace=True)
        return candidates_df
    else:
        return pd.DataFrame()


def create_pileup(in_dct):
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}

    if in_dct['mode']=='testing':
        test_sites_df=get_candidates(in_dct)
        if test_sites_df.shape[0]==0:
            return None

        test_sites_df.set_index('POS',drop=False,inplace=True)
        d={}
        for v_pos in test_sites_df.POS:

            res=create_pileup_image(in_dct,v_pos=v_pos)
            if res:
                d[v_pos]=(test_sites_df['REF'].loc[v_pos],res[0])

        return d

    elif in_dct['mode']=='common':
        d={}
        fastafile=pysam.FastaFile(in_dct['fasta_path'])
        for v_pos in in_dct['list']:
            r=fastafile.fetch(in_dct['chrom'],v_pos-1,v_pos)
            if r in ['A','G','T','C']:
                res=create_pileup_image(in_dct,v_pos=v_pos)
                if res:
                    d[v_pos]=(mapping[r],res[0])

        return d
        
        
        
def create_pileup_image(dct,v_pos=None):
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
        
        d={x:{} for x in range(v_start,v_end+1)}
        features={x:{} for x in range(v_start,v_end)}
        try:
            rlist=[mapping[s] for s in fastafile.fetch(fasta_chr,v_start-1,v_end)]
            
        except KeyError:
            #print('KeyError at %d' %v_pos, flush=True)
            return None
        
        for pcol in samfile.pileup(chrom,v_pos-1,v_pos,min_base_quality=0,stepper='nofilter',\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            if pcol.get_num_aligned()<4:
                return None
            
        for pcol in samfile.pileup(chrom,v_start-1,v_end,min_base_quality=0,stepper='nofilter',\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
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
        tmp=np.dstack([(p_mat==i) for i in range(4)])
        
        #data=np.multiply(tmp,tmp2).astype(np.int8)
        try:
            ref_match=2*ref_match-1
            f_mat=np.array(f_df)*ref_match
            data=np.dstack((tmp,f_mat))
        except ValueError:
            #print('ValueError at %d' %v_pos, flush=True)
            return None
        
        if data.shape[0]<du_lim:
            tmp=np.zeros((du_lim-data.shape[0],data.shape[1],data.shape[2]))
            data=np.vstack((data,tmp))
            rnames+=['Buff' for i in range(tmp.shape[0])]
            
        return (data.astype(np.int8),rnames)
    except AssertionError:
        #print('AssertionError at %d' %v_pos, flush=True)
        return None
        
def generate(params,mode='training'):
    cores=params['cpu']
    mode=params['mode']
    if mode=='training':
        #print('starting pileups',flush=True)
        
        pool = mp.Pool(processes=cores)
        fname='%s_pileups' %params['chrom']
        true_fname=os.path.join(params['out_path'],fname+'_pos')
        p_file=open(true_fname , "w")
        
        neg_file_list={}
        for i in [0,5,10,15,20,25]:
            tmp_name=os.path.join(params['out_path'],'%s_neg_%d' %(fname,i))
            neg_file_list[i]=open(tmp_name , "w")
        
        start,end=params['start'],params['end']
        
        
        
        for mbase in range(start,end,int(1e6)):
            print('starting pool:'+str(mbase),flush=True)
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
                neg_pile_list=result[1]
                for i in [0,5,10,15,20,25]:
                    neg_pile=neg_pile_list[i]
                    if neg_pile:
                        for data in neg_pile:
                                pos,gt,allele,ref,mat=data
                                ft=mat[:,:,-1].reshape(-1)
                                mat=mat[:,:,:-1].reshape(-1)
                                s='%s%d%d%d%d%s%s' %((11-len(str(pos)))*'0',pos,gt,allele,ref, ''.join(mat.astype('<U1')),''.join([(3-len(x))*' '+x for x in ft.astype('<U3')]))
                                neg_file_list[i].write(s)

                if result[0]:
                    for data in result[0]:
                        pos,gt,allele,ref,mat=data

                        ft=mat[:,:,-1].reshape(-1)
                        mat=mat[:,:,:-1].reshape(-1)
                        s='%s%d%d%d%d%s%s' %((11-len(str(pos)))*'0',pos,gt,allele,ref, ''.join(mat.astype('<U1')),''.join([(3-len(x))*' '+x for x in ft.astype('<U3')]))
                        p_file.write(s)

            
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)
            
    elif mode=='testing':
        print('starting pileups',flush=True)
        pool = mp.Pool(processes=cores)
        fname='%s_pileups_test' %(params['chrom'])
        fname=os.path.join(params['out_path'],fname)
        file=open(fname , "w")
        start,end=params['start'],params['end']
        for mbase in range(start,end,int(1e4)):
            print('starting pool:'+str(mbase),flush=True)
            t=time.time()                

            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(1e4)),1000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=k+1000
                in_dict_list.append(d)
            results = pool.map(create_pileup, in_dict_list)
            for res_dict in results:
                if res_dict:
                    for k,p in res_dict.items():
                        ref,mat=p

                        ft=mat[:,:,-1].reshape(-1)
                        mat=mat[:,:,:-1].reshape(-1)
                        s='%s%d%d%s%s' %((11-len(str(k)))*'0',k,ref, ''.join(mat.astype('<U1')),''.join([(3-len(x))*' '+x for x in ft.astype('<U3')]))
                        file.write(s)

            
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)
    
    elif mode=='common':
        with open(params['vcf_path'],'r') as file:
            content=[int(x.rstrip('\n')) for x in file]
        
            
        ind_list=[content[i*(len(content)//10):(i+1)*(len(content)//10)] for i in range(10+(len(content)%10!=0))]
        pool = mp.Pool(processes=cores)
        print('starting pileups',flush=True)
        for c,ind in enumerate(ind_list):
            t=time.time()
            fname='%s_known_sites_pileups_%d' %(params['chrom'],c+1)
            print('starting: %s' %fname,flush=True)
            fname=os.path.join(params['out_path'],fname)
            
            pool_list=[ind[i*(len(ind)//cores):(i+1)*(len(ind)//cores)] for i in range(cores+(len(content)%cores!=0))]
            
            in_dict_list=[]
            for k in pool_list:
                d = copy.deepcopy(params)    
                d['list']=k
                in_dict_list.append(d)
            results = pool.map(create_pileup, in_dict_list)
            
                
            with open(fname,'w') as file:
                for res_dict in results:
                    if res_dict:
                        for k,p in res_dict.items():
                            ref,mat=p

                            ft=mat[:,:,-1].reshape(-1)
                            mat=mat[:,:,:-1].reshape(-1)
                            s='%s%d%d%s%s' %((11-len(str(k)))*'0',k,ref, ''.join(mat.astype('<U1')),''.join([(3-len(x))*' '+x for x in ft.astype('<U3')]))
                            file.write(s)
            
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool',flush=True)
            

        
        


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
    
    args = parser.parse_args()
    
    if len(args.region.split(':'))==2:
        chrom,region=args.region.split(':')
        start,end=int(region.split('-')[0]),int(region.split('-')[1])
        
    else:
        chrom=args.region.split(':')[0]
        start,end=1,chrom_length[chrom]
        
    in_dict={'mode':args.mode, 'chrom':chrom,'start':start,'end':end,\
         'sam_path':args.bam,'fasta_path':args.ref,'vcf_path':args.vcf,'out_path':args.output,'window':args.window,'depth':args.depth,'threshold':args.threshold,'cpu':args.cpu}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
    
