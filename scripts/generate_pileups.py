import sys,pysam, time,os,re,copy,argparse,gzip,itertools
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree


def get_training_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    vcf_path=dct['vcf_path']
    bed_path=dct['bed']
    window=dct['window']
    nbr_size=dct['nbr_size']
    cnd_df=dct['cand_list']
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
    
    tr_output=[]
    output={'pos':[],0:[],5:[],10:[],15:[],20:[],25:[]}
    
    if bed_path:
        bed_file=os.path.join(bed_path,'%s.bed' %chrom)
        with open(bed_file) as file:
            content=[x.rstrip('\n') for x in file]

        content=[x.split('\t')[1:] for x in content]
        content=[(int(x[0]),int(x[1])) for x in content]
        t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)


    samfile = pysam.Samfile(sam_path, "rb")
    
    bcf_in = VariantFile(vcf_path)  # auto-detect input format
    fastafile=pysam.FastaFile(fasta_path)
    
    gt_map={(0,0):0, (1,1):0, (2,2):0, (1,2):1, (2,1):1, (0,1):1, (1,0):1, (0,2):1,(2,0):1}
    tr_pos={}
    for rec in bcf_in.fetch(chrom,start,end+1):
        tr_pos[rec.pos]=(gt_map[rec.samples.items()[0][1].get('GT')],mapping[rec.alleles[1]])
        

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
    
    pileup_dict={}
    
    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            
            if ref_df.loc[pcol.pos+1].ref!='*':
                    
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()
                n=pcol.get_num_aligned()

                pileup_dict[pcol.pos+1]={n:s for (n,s) in zip(name,seq)}


                if bed_path:
                    if not t[pcol.pos+1]:
                        continue

                if pcol.pos+1 in tr_pos.keys() and n>=dct['mincov'] :
                    output['pos'].append(pcol.pos+1)
                    
                if n>=dct['mincov'] and pcol.pos+1>=start and pcol.pos+1<=end and pcol.pos+1 not in tr_pos:
                    alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=ref_df.loc[pcol.pos+1].ref and x[0] in 'AGTC')]+[0])/n

                    if alt_freq<0.05:
                        output[0].append(pcol.pos+1)

                    elif 0.05<=alt_freq<0.10:
                        output[5].append(pcol.pos+1)
                    elif 0.10<=alt_freq<0.15:
                        output[10].append(pcol.pos+1)
                    elif 0.15<=alt_freq<0.20:
                        output[15].append(pcol.pos+1)
                    elif 0.20<=alt_freq<0.25:
                        output[20].append(pcol.pos+1)
                    elif 0.25<=alt_freq:
                        output[25].append(pcol.pos+1)
                    
    
    pileup_list={'pos':[],0:[],5:[],10:[],15:[],20:[],25:[]}
    
    ref_df[['ref']]=ref_df[['ref']].applymap(lambda x:mapping[x])
    
    pileup_df=pd.DataFrame.from_dict(pileup_dict)
    pileup_df.dropna(axis=1, how='all',inplace=True)

    np.random.seed(76)
    
    if output['pos']:
        tr_len=len(output['pos'])    
    else:
        tr_len=1e16
        
    sizes={0:2*tr_len, 5:tr_len//3,10:tr_len//3,15:tr_len//3, 20:tr_len//2, 25:tr_len}

    for i in ['pos',0,5,10,15,20,25]:
        pos_list=output[i]
        
        if pos_list:
            if i!='pos':
                if sizes[i]<len(output[i]):
                    perm=np.random.permutation(sizes[i])
                    pos_list=np.take(pos_list,perm,axis=0)
        
            for v_pos in pos_list:
                p_df=pileup_df.reindex(list(range(v_pos-window,v_pos+window+1)),axis=1)
                ls=cnd_df[(cnd_df['pos']<v_pos+20000) & (cnd_df['pos']>v_pos-20000)].pos
                #p_df.dropna(axis=1, how='all',inplace=True)

                ls1=[p for p in ls if p<v_pos-window][-nbr_size:]
                ls2=[p for p in ls if p>v_pos+window][:nbr_size]

                nbr1_dict={}
                nbr2_dict={}

                rlist1=[]
                rlist2=[]
                
                for nb_pos in ls1:
                            nbr1_dict[nb_pos]={n:s for (n,s) in zip(cnd_df.loc[nb_pos].names.split(':'),cnd_df.loc[nb_pos].seq) if n in p_df.index}
                            rlist1.append(cnd_df.loc[nb_pos].ref)

                for nb_pos in ls2:
                            nbr2_dict[nb_pos]={n:s for (n,s) in zip(cnd_df.loc[nb_pos].names.split(':'),cnd_df.loc[nb_pos].seq) if n in p_df.index}
                            rlist2.append(cnd_df.loc[nb_pos].ref)

                
                total_rlist=rlist1+ref_df.loc[p_df.columns.to_list()].ref.to_list()+rlist2
                
                nbr1_df=pd.DataFrame.from_dict(nbr1_dict)
                nbr2_df=pd.DataFrame.from_dict(nbr2_dict)

                p_df=pd.concat([nbr1_df,p_df,nbr2_df],axis=1,sort=False)

                p_df.dropna(subset=[v_pos],inplace=True)

                p_df['counter']=p_df[v_pos]
                new_df=p_df.groupby('counter').agg(['value_counts']).reindex(list(itertools.product('AGTC','AGTC'))).fillna(0)
                mat=np.array(new_df).reshape([4,4,new_df.shape[1]]).swapaxes(1,2)

                total_ref=np.eye(5)[total_rlist]
                total_ref=total_ref[np.newaxis,:]
                
                total_ref[:,4]=0
                
                mat=np.dstack([mat,np.zeros([4,mat.shape[1]])+np.eye(4)[ref_df.loc[v_pos].ref][:,np.newaxis]])
                
                data=np.vstack([total_ref,np.multiply(mat,1-2*total_ref)])
                data=np.hstack([np.zeros([5,nbr_size-len(ls1),5]),data,np.zeros([5,nbr_size-len(ls2),5])]).astype(np.int8)

                if i=='pos':
                    pileup_list[i].append((v_pos,tr_pos[v_pos][0],tr_pos[v_pos][1],ref_df.loc[v_pos].ref,data))
                else:
                    pileup_list[i].append((v_pos,0,ref_df.loc[v_pos].ref,ref_df.loc[v_pos].ref,data))
               

    return pileup_list




def get_testing_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    threshold=dct['threshold']
    bed_path=dct['bed']
    window=dct['window']
    nbr_size=dct['nbr_size']
    cnd_df=dct['cand_list']
    
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
    
    if bed_path:
        bed_file=os.path.join(bed_path,'%s.bed' %chrom)
        with open(bed_file) as file:
            content=[x.rstrip('\n') for x in file]

        content=[x.split('\t')[1:] for x in content]
        content=[(int(x[0]),int(x[1])) for x in content]
        t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)


    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
    
    pileup_dict={}
    
    output={5:[],10:[],15:[],20:[]}
    
    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            if ref_df.loc[pcol.pos+1].ref!='*':
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()
                n=pcol.get_num_aligned()

                pileup_dict[pcol.pos+1]={n:s for (n,s) in zip(name,seq)}


                if bed_path:
                    if not t[pcol.pos+1]:
                        continue


                if n>=dct['mincov'] and pcol.pos+1>=start and pcol.pos+1<=end:
                    alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=ref_df.loc[pcol.pos+1].ref and x[0] in 'AGTC')]+[0])/n

                    if 0.5<=alt_freq<0.10 and 5 in threshold:
                            output[5].append(pcol.pos+1)
                    
                    elif 0.10<=alt_freq<0.15 and 10 in threshold:
                            output[10].append(pcol.pos+1)
                        
                    elif 0.15<=alt_freq<0.20 and 15 in threshold:
                            output[15].append(pcol.pos+1)
                            
                    elif 0.20<=alt_freq and 20 in threshold:
                            output[20].append(pcol.pos+1)
                    
    pileup_list={x:[] for x in threshold}
    
    ref_df[['ref']]=ref_df[['ref']].applymap(lambda x:mapping[x])
    
    pileup_df=pd.DataFrame.from_dict(pileup_dict)
    pileup_df.dropna(axis=1, how='all',inplace=True)
    
    for i in threshold:
        pos_list=output[i]
        
        if pos_list:

            for v_pos in pos_list:
                p_df=pileup_df.reindex(list(range(v_pos-window,v_pos+window+1)),axis=1)
                ls=cnd_df[(cnd_df['pos']<v_pos+20000) & (cnd_df['pos']>v_pos-20000)].pos
                #p_df.dropna(axis=1, how='all',inplace=True)

                ls1=[p for p in ls if p<v_pos-window][-nbr_size:]
                ls2=[p for p in ls if p>v_pos+window][:nbr_size]


                nbr1_dict={}
                nbr2_dict={}

                rlist1=[]
                rlist2=[]
                for nb_pos in ls1:
                            nbr1_dict[nb_pos]={n:s for (n,s) in zip(cnd_df.loc[nb_pos].names.split(':'),cnd_df.loc[nb_pos].seq) if n in p_df.index}
                            rlist1.append(cnd_df.loc[nb_pos].ref)

                for nb_pos in ls2:
                            nbr2_dict[nb_pos]={n:s for (n,s) in zip(cnd_df.loc[nb_pos].names.split(':'),cnd_df.loc[nb_pos].seq) if n in p_df.index}
                            rlist2.append(cnd_df.loc[nb_pos].ref)

                total_rlist=rlist1+ref_df.loc[p_df.columns.to_list()].ref.to_list()+rlist2
                nbr1_df=pd.DataFrame.from_dict(nbr1_dict)
                nbr2_df=pd.DataFrame.from_dict(nbr2_dict)

                p_df=pd.concat([nbr1_df,p_df,nbr2_df],axis=1,sort=False)

                p_df.dropna(subset=[v_pos],inplace=True)

                p_df['counter']=p_df[v_pos]
                
                new_df=p_df.groupby('counter').agg(['value_counts']).reindex(list(itertools.product('AGTC','AGTC'))).fillna(0)
                mat=np.array(new_df).reshape([4,4,new_df.shape[1]]).swapaxes(1,2)

                total_ref=np.eye(5)[total_rlist]
                total_ref=total_ref[np.newaxis,:]
                
                total_ref[:,4]=0
                
                mat=np.dstack([mat,np.zeros([4,mat.shape[1]])+np.eye(4)[ref_df.loc[v_pos].ref][:,np.newaxis]])
                
                data=np.vstack([total_ref,np.multiply(mat,1-2*total_ref)])
                data=np.hstack([np.zeros([5,nbr_size-len(ls1),5]),data,np.zeros([5,nbr_size-len(ls2),5])]).astype(np.int8)
                
                pileup_list[i].append((v_pos,ref_df.loc[v_pos].ref,data))


    return pileup_list
    
    
def redo_candidates(dct):
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    threshold=dct['threshold']
    bed_path=dct['bed']
    window=dct['window']
    nbr_size=dct['nbr_size']
    cnd_df=dct['cand_list']

    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
    
    pileup_dict={}
    
    output=[]
    
    bcf_in = VariantFile(dct['vcf_path'])
    for rec in bcf_in.fetch(chrom,start,end):
        output.append(rec.pos)
                
    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            if pcol.pos+1 in ref_df.pos:
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()
                n=pcol.get_num_aligned()

                pileup_dict[pcol.pos+1]={n:s for (n,s) in zip(name,seq)}

    if output:
        ts_list=[]
        mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
        ref_df[['ref']]=ref_df[['ref']].applymap(lambda x:mapping[x])
        pileup_df=pd.DataFrame.from_dict(pileup_dict)
        pileup_df.dropna(axis=1, how='all',inplace=True)

        for v_pos in output:
            p_df=pileup_df.reindex(list(range(v_pos-window,v_pos+window+1)),axis=1)
            #p_df.dropna(axis=1, how='all',inplace=True)

            ls_10=cnd_df[(cnd_df['gt']=='1/0') & (cnd_df['pos']<v_pos+20000) & (cnd_df['pos']>v_pos-20000)].pos
            ls_01=cnd_df[(cnd_df['gt']=='0/1') & (cnd_df['pos']<v_pos+20000) & (cnd_df['pos']>v_pos-20000)].pos
            
            
            ls1=[p for p in ls_10 if p<v_pos-window][-nbr_size:]
            
            if len(ls1)<nbr_size:
                    ls1=ls1+[p for p in ls_01 if p<v_pos-window][-(nbr_size-len(ls1)):]
                    ls1.sort()
                    
            
            ls2=[p for p in ls_10 if p>v_pos+window][:nbr_size]
            
            if len(ls2)<nbr_size:
                    ls2=ls2+[p for p in ls_01 if p>v_pos+window][-(nbr_size-len(ls2)):]
                    ls2.sort()
                    
            nbr1_dict={}
            nbr2_dict={}
            
            rlist1=[]
            rlist2=[]
            
            for nb_pos in ls1:
                        nbr1_dict[nb_pos]={n:s for (n,s) in zip(cnd_df.loc[nb_pos].names.split(':'), cnd_df.loc[nb_pos].seq) if n in p_df.index}
                        rlist1.append(cnd_df.loc[nb_pos].ref)
            

                
            for nb_pos in ls2:
                        nbr2_dict[nb_pos]={n:s for (n,s) in zip(cnd_df.loc[nb_pos].names.split(':'),cnd_df.loc[nb_pos].seq) if n in p_df.index}
                        rlist2.append(cnd_df.loc[nb_pos].ref)
            
            total_rlist=rlist1+ref_df.loc[p_df.columns.to_list()].ref.to_list()+rlist2
            nbr1_df=pd.DataFrame.from_dict(nbr1_dict)
            nbr2_df=pd.DataFrame.from_dict(nbr2_dict)

            
            p_df=pd.concat([nbr1_df,p_df,nbr2_df],axis=1,sort=False)

            p_df.dropna(subset=[v_pos],inplace=True)
            
            p_df['counter']=p_df[v_pos]
            new_df=p_df.groupby('counter').agg(['value_counts']).reindex(list(itertools.product('AGTC','AGTC'))).fillna(0)
            mat=np.array(new_df).reshape([4,4,new_df.shape[1]]).swapaxes(1,2)

            total_ref=np.eye(5)[total_rlist]
            total_ref=total_ref[np.newaxis,:]

            total_ref[:,4]=0

            mat=np.dstack([mat,np.zeros([4,mat.shape[1]])+np.eye(4)[ref_df.loc[v_pos].ref][:,np.newaxis]])

            data=np.vstack([total_ref,np.multiply(mat,1-2*total_ref)])
            data=np.hstack([np.zeros([5,nbr_size-len(ls1),5]),data,np.zeros([5,nbr_size-len(ls2),5])]).astype(np.int8)

            ts_list.append((v_pos,ref_df.loc[v_pos].ref,data))

            
        p_df,pileup_df,ref_df,pileup_dict=None,None,None,None
        return ts_list
    else:
        return []

def generate(params,mode='training'):
    print('adfasd',flush=True)
    cores=params['cpu']
    mode=params['mode']
    chrom=params['chrom']
    threshold=params['threshold']
    start,end=params['start'],params['end']
    
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
    
    if mode in ['training','train']: #'pos,ref,seq,names'
        print('starting training pileups',flush=True)
        cnd_df=pd.read_csv(os.path.join(params['cnd_path'],'%s.pileups.neighbors.train' %chrom))
        cnd_df.rename(columns={0:'pos',1:'depth',2:'freq',3:'ref',4:'seq',5:'names'},inplace=True)
        cnd_df=cnd_df[(cnd_df['pos']>=start-20000)& (cnd_df['pos']<=end+20000)]
        cnd_df.set_index('pos',inplace=True,drop=False)
        cnd_df[['ref']]=cnd_df[['ref']].applymap(lambda x:mapping[x])
        params['cand_list']=cnd_df


        pool = mp.Pool(processes=cores)
        fname='%s.pileups.' %params['chrom']
        file_list={}
        for i in ['pos',0,5,10,15,20,25]:
            suffix='pos' if i == 'pos' else 'neg.%d' %i
            tmp_name=os.path.join(params['out_path'],'%s.%s' %(fname,suffix)) 
            file_list[i]=open(tmp_name , "w")
        
        for mbase in range(start,end,int(1e7)):
            print('starting pool:'+str(mbase),flush=True)
            t=time.time()                

            in_dict_list=[]

            for k in range(mbase,min(end,mbase+int(1e7)),100000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(k+100000,end)
                in_dict_list.append(d)
            results_dict = pool.map(get_training_candidates, in_dict_list)
                
            for result in results_dict:
                for i in ['pos',0,5,10,15,20,25]:
                    pileups=result[i]
                    if pileups:
                        for data in pileups:
                                pos,gt,allele,ref,mat=data
                                mat=mat.reshape(-1)
                                s='%s%d%d%d%d%s' %((11-len(str(pos)))*'0',pos,gt,allele,ref,''.join([(6-len(x))*' '+x for x in mat.astype('<U6')]))
                                file_list[i].write(s)
            
            results_dict=None
            elapsed=time.time()-t
            
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)                    

                
    elif mode in ['testing','test']:#'pos,depth,freq,ref,seq,names'
        cnd_df=pd.read_csv(os.path.join(params['cnd_path'],'%s.pileups.neighbors.test' %chrom))
        
        cnd_df=cnd_df[(cnd_df['pos']>=start-20000)& (cnd_df['pos']<=end+20000)]

        cnd_df=cnd_df[cnd_df['ref'].map(lambda x: x in ['A','G','T','C'])]

        cnd_df[['ref']]=cnd_df[['ref']].applymap(lambda x:mapping[x])


        cnd_df.set_index('pos',inplace=True,drop=False)
        params['cand_list']=cnd_df
        
        print('starting pileups',flush=True)
        pool = mp.Pool(processes=cores)
        fname='%s.pileups.' %params['chrom']
        file_list={}
        for i in threshold:
            suffix='test.%d' %i
            tmp_name=os.path.join(params['out_path'],'%s%s' %(fname,suffix)) 
            file_list[i]=open(tmp_name , "w")

        start,end=params['start'],params['end']
        
        for mbase in range(start,end,int(1e7)):
            print('starting pool:'+str(mbase),flush=True)
            t=time.time()                

            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(1e7)),100000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(end,k+100000)
                in_dict_list.append(d)
            results_dict = pool.map(get_testing_candidates, in_dict_list)
            
            for result in results_dict:
                for i in threshold:
                    pileups=result[i]
                    if pileups:
                        for data in pileups:
                            
                            pos,ref,mat=data
                            mat=mat.reshape(-1)
                            s='%s%d%d%s' %((11-len(str(pos)))*'0',pos,ref,''.join([(6-len(x))*' '+x for x in mat.astype('<U6')]))
                            file_list[i].write(s)

            results_dict=None
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)
            
            
    elif mode=='redo':#'pos,ref,gt,seq,names'
        cnd_df=pd.read_csv(os.path.join(params['cnd_path'],'%s.pileups.neighbors.redo' %chrom))
        cnd_df=cnd_df[(cnd_df['pos']>=start-20000)& (cnd_df['pos']<=end+20000)]

        cnd_df=cnd_df[cnd_df['ref'].map(lambda x: x in ['A','G','T','C'])]

        
        cnd_df[['ref']]=cnd_df[['ref']].applymap(lambda x:mapping[x])

        cnd_df.set_index('pos',inplace=True,drop=False)
        

        params['cand_list']=cnd_df
        
        print('starting pileups',flush=True)
        pool = mp.Pool(processes=cores)
        fname='%s.pileups.test.redov3' %params['chrom']
        fname=os.path.join(params['out_path'],fname)
        file=open(fname , "w")
        start,end=params['start'],params['end']
        for mbase in range(start,end,int(1e7)):
            print('starting pool:'+str(mbase),flush=True)
            t=time.time()                

            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(1e7)),100000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(end,k+100000)
                in_dict_list.append(d)
            results_dict = pool.map(redo_candidates, in_dict_list)
            
            for result in results_dict:
                if result:
                    for data in result:
                        pos,ref,mat=data
                        mat=mat.reshape(-1)
                        s='%s%d%d%s' %((11-len(str(pos)))*'0',pos,ref,''.join([(6-len(x))*' '+x for x in mat.astype('<U6')]))
                        file.write(s)

            results_dict=None
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)


if __name__ == '__main__':
    chrom_length={'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555, 'chr5':181538259, 'chr6':170805979, \
             'chr7':159345973, 'chr8':145138636, 'chr9':138394717, 'chr10':133797422, 'chr11':135086622, 'chr12':133275309,\
             'chr13':114364328, 'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,\
             'chr19':58617616, 'chr20':64444167, 'chr21':46709983, 'chr22':50818468, 'chrX':156040895, 'chrY':57227415}
    parser = argparse.ArgumentParser()

    #-r chromosome region   -m mode   -bam bam file   -ref reference file   -vcf ground truth variants   -o output path
    parser.add_argument("-chrom", "--chrom", help="Chromosome region")
    parser.add_argument("-m", "--mode", help="Mode")
    parser.add_argument("-bam", "--bam", help="Bam file")
    parser.add_argument("-ref", "--ref", help="Size")
    parser.add_argument("-vcf", "--vcf", help="Ground truth variants")
    parser.add_argument("-o", "--output", help="Output path")
    parser.add_argument("-w", "--window", help="Window",type=int)
    parser.add_argument("-cpu", "--cpu", help="CPUs",type=int)
    parser.add_argument("-d", "--depth", help="Depth",type=int)
    parser.add_argument("-t", "--threshold", help="Threshold")
    parser.add_argument("-bed", "--bed", help="BED file")
    parser.add_argument("-cnd", "--cnd_path", help="Candidates file")
    parser.add_argument("-mincov", "--mincov", help="min coverage",type=int)
    parser.add_argument("-start", "--start", help="start",type=int)
    parser.add_argument("-end", "--end", help="end",type=int)
    parser.add_argument("-nbr", "--nbr_size", help="Number of neighboring het SNPs",type=int)

    
    args = parser.parse_args()
    
    chrom=args.chrom
    
    if not args.end:
        end=chrom_length[chrom]
    else:
        end=args.end
    
    if not args.start:
        start=1
    else:
        start=args.start
    if not args.threshold:
        threshold=[15,20]
    else:
        threshold=[int(x) for x in args.threshold.split(':')]
    in_dict={'mode':args.mode, 'chrom':chrom,'start':start,'end':end,\
         'sam_path':args.bam, 'fasta_path':args.ref, 'vcf_path':args.vcf,\
             'out_path':args.output, 'window':args.window, 'depth':args.depth,\
             'threshold':threshold, 'cpu':args.cpu, 'bed':args.bed,'cnd_path':args.cnd_path,\
            'mincov':args.mincov, 'nbr_size':args.nbr_size}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
    

    
    
