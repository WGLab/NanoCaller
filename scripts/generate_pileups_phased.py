import sys,pysam, time,os,re,copy,argparse,gzip,itertools
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree


def get_training_candidates(dct):
    du_lim=32
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    vcf_path=dct['vcf_path']
    bed_path=dct['bed']
    window=dct['window']

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
    
    fastafile=pysam.FastaFile(fasta_path)
    
    tr_df=dct['tr_df'][0]
    total_pos=dct['tr_df'][1]
    tr_df=tr_df[(tr_df.pos>=start)&(tr_df.pos<=end)]
    
    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
    
    hap_df=pd.read_csv('%s.%s.bam.info' %(dct['hap_path'],chrom), delim_whitespace=True)

    hap_reads_0=set(hap_df[(hap_df.haplotype==0)]['#readname'])
    hap_reads_1=set(hap_df[(hap_df.haplotype==1)]['#readname'])

    flag_dict={}
    #q_dict={}
    #a_score={}
    for pread in samfile.fetch(chrom,max(0,start-1-30),end+30):
        flag_dict[pread.qname]=(pread.flag &0x10)//16

    pileup_dict={}
    features={}
    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            
            if ref_df.loc[pcol.pos+1].ref in 'AGTC':
                    
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)]).upper()
                
                name=pcol.get_query_names()
                qual=pcol.get_query_qualities()
                dp=len(seq)
                n=len(seq.replace('*',''))
                
                pileup_dict[pcol.pos+1]={n:s for (n,s) in zip(name,seq)}
                features[pcol.pos+1]={n:s for (n,s) in zip(name,qual)}
                
                if bed_path:
                    if not t[pcol.pos+1]:
                        continue

                if pcol.pos+1 in tr_df.pos and n>=dct['mincov'] :
                    output['pos'].append(pcol.pos+1)
                    
                elif pcol.pos+1 not in total_pos and np.random.randint(4)==1:
                    if dp>=dct['mincov'] and pcol.pos+1>=start and pcol.pos+1<=end and pcol.pos+1 not in tr_df.pos:
                        alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=ref_df.loc[pcol.pos+1].ref and x[0] in 'AGTC')]+[0])/n if n!=0 else 0

                        if 0<alt_freq<0.15:
                            output[0].append(pcol.pos+1)

                        elif 0.15<=alt_freq:
                            output[15].append(pcol.pos+1)
                    
    
    pileup_list={'pos':[],0:[],15:[]}
    
    ref_df[['ref']]=ref_df[['ref']].applymap(lambda x:mapping[x])

    np.random.seed(76)
    
    if output['pos']:
        tr_len=len(output['pos'])    
    else:
        tr_len=1e16
        
    sizes={0:tr_len//4, 15:tr_len//4}

    for i in ['pos',0,15]:
        pos_list=output[i]
        
        if pos_list:
            if i!='pos':
                if sizes[i]<len(output[i]):
                    perm=np.random.permutation(sizes[i])
                    pos_list=np.take(pos_list,perm,axis=0)
        
            for v_pos in pos_list:
                
                rlist=np.array([ref_df.loc[x].ref for x in list(range(v_pos-window,v_pos+window+1)) if x in pileup_dict.keys()])

                p_df=pd.DataFrame({x:pileup_dict[x] for x in list(range(v_pos-window,v_pos+window+1)) if x in pileup_dict.keys()})
                f_df=pd.DataFrame({x:features[x] for x in list(range(v_pos-window,v_pos+window+1)) if x in features.keys()})

                p_df.fillna('N',inplace=True)
                f_df.fillna(0,inplace=True)
                
                p_df=p_df.applymap(lambda x: mapping[x])

                p_df_0=p_df.reindex(hap_reads_0&set(p_df.index))
                if p_df_0.shape[0]>du_lim:
                        p_df_0=p_df_0.sample(n=du_lim,replace=False, random_state=1)

                p_df_0=p_df_0.assign(f = p_df_0.index.map(lambda x:flag_dict[x]))
                p_df_0.sort_values([v_pos,'f'],inplace=True,ascending=False)
                p_df_0_flag=p_df_0.f
                p_df_0.drop('f', axis=1,inplace=True)
                p_df_0.dropna(subset=[v_pos],inplace=True)
                f_df_0=f_df.reindex(p_df_0.index)
                f_mat_0=np.array(f_df_0)
                p_mat_0=np.array(p_df_0)
                ref_0=(p_mat_0==rlist)
                p_mat_0=np.eye(6)[p_mat_0]
                p_mat_0=p_mat_0[:,:,:4]
                strand_0=np.array(p_df_0_flag)[:,np.newaxis]+np.zeros(p_mat_0.shape[:-1])
                data_0=np.dstack([p_mat_0,ref_0,f_mat_0,strand_0])

                p_df_1=p_df.reindex(hap_reads_1&set(p_df.index))
                if p_df_1.shape[0]>du_lim:
                        p_df_1=p_df_1.sample(n=du_lim,replace=False, random_state=1)
                p_df_1=p_df_1.assign(f = p_df_1.index.map(lambda x:flag_dict[x]))
                p_df_1.sort_values([v_pos,'f'],inplace=True,ascending=False)
                p_df_1_flag=p_df_1.f
                p_df_1.drop('f', axis=1,inplace=True)
                p_df_1.dropna(subset=[v_pos],inplace=True)
                f_df_1=f_df.reindex(p_df_1.index)
                f_mat_1=np.array(f_df_1)
                p_mat_1=np.array(p_df_1)
                ref_1=(p_mat_1==rlist)
                p_mat_1=np.eye(6)[p_mat_1]
                p_mat_1=p_mat_1[:,:,:4]
                strand_1=np.array(p_df_1_flag)[:,np.newaxis]+np.zeros(p_mat_1.shape[:-1])
                data_1=np.dstack([p_mat_1,ref_1,f_mat_1,strand_1])

                tmp0=True
                if data_0.shape[0]<2:
                    tmp0=False
                    
                elif data_0.shape[0]<du_lim:
                        tmp=np.zeros((du_lim-data_0.shape[0],data_0.shape[1],data_0.shape[2]))
                        data_0=np.vstack((data_0,tmp))

                tmp1=True
                if data_1.shape[0]<2:
                    tmp1=False
                elif data_1.shape[0]<du_lim:
                        tmp=np.zeros((du_lim-data_1.shape[0],data_1.shape[1],data_1.shape[2]))
                        data_1=np.vstack((data_1,tmp))


                if i=='pos':
                    if tmp0:
                        pileup_list[i].append((v_pos,tr_df.loc[v_pos][0],ref_df.loc[v_pos].ref, data_0.astype(np.int16)))
                    
                    if tmp1:
                        pileup_list[i].append((v_pos,tr_df.loc[v_pos][1],ref_df.loc[v_pos].ref, data_1.astype(np.int16)))
                else:
                    if tmp0:
                        pileup_list[i].append((v_pos,ref_df.loc[v_pos].ref,ref_df.loc[v_pos].ref, data_0.astype(np.int16)))
                    
                    if tmp1:
                        pileup_list[i].append((v_pos,ref_df.loc[v_pos].ref, ref_df.loc[v_pos].ref, data_1.astype(np.int16)))

    return pileup_list

def get_testing_candidates(dct):
    du_lim=32
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    vcf_path=dct['vcf_path']
    bed_path=dct['bed']
    window=dct['window']

    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
    
    output={'pos':[],0:[],5:[],10:[],15:[],20:[],25:[]}
    
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
    
    hap_df=pd.read_csv('%s.%s.bam.info' %(dct['hap_path'],chrom), delim_whitespace=True)

    hap_reads_0=set(hap_df[(hap_df.haplotype==0)]['#readname'])
    hap_reads_1=set(hap_df[(hap_df.haplotype==1)]['#readname'])

    flag_dict={}
    #q_dict={}
    #a_score={}
    for pread in samfile.fetch(chrom,max(0,start-1-30),end+30):
        flag_dict[pread.qname]=(pread.flag &0x10)//16

    pileup_dict={'hap0':{},'hap1':{}}
    features={'hap0':{},'hap1':{}}
    
    output={'hap0':{},'hap1':{}}
    
    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            
            if ref_df.loc[pcol.pos+1].ref in 'AGTC':
                r=ref_df.loc[pcol.pos+1].ref
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)]).upper()
                
                name=pcol.get_query_names()
                qual=pcol.get_query_qualities()
                
                pileup_dict['hap0'][pcol.pos+1]={n:s for (n,s) in zip(name,seq) if n in hap_reads_0}
                features['hap0'][pcol.pos+1]={n:s for (n,s) in zip(name,qual) if n in hap_reads_0}
                
                pileup_dict['hap1'][pcol.pos+1]={n:s for (n,s) in zip(name,seq) if n in hap_reads_1}
                features['hap1'][pcol.pos+1]={n:s for (n,s) in zip(name,qual) if n in hap_reads_1}
                
                if bed_path:
                    if not t[pcol.pos+1]:
                        continue
                        
                hap0_bases=''.join(pileup_dict['hap0'][pcol.pos+1]).replace('*','')
                hap1_bases=''.join(pileup_dict['hap1'][pcol.pos+1]).replace('*','')
                seq=seq.replace('*','')

                tot_freq=max([x[1] for x in Counter(seq).items() if x[0]!=r]+[0])/len(seq) if len(seq)>0 else 0

                hap0_freq=max([x[1] for x in Counter(hap0_bases).items() if x[0]!=r]+[0])/len(hap0_bases) if len(hap0_bases)>0 else 0
                hap1_freq=max([x[1] for x in Counter(hap1_bases).items() if x[0]!=r]+[0])/len(hap1_bases) if len(hap1_bases)>0 else 0

                if hap0_freq>=0.4 and len(pileup_dict['hap0'][pcol.pos+1])>=dct['mincov']:
                    output['hap0'][pcol.pos+1]=(len(pileup_dict['hap0'][pcol.pos+1]),hap0_freq)
                if hap1_freq>=0.4 and len(pileup_dict['hap1'][pcol.pos+1])>=dct['mincov']:
                    output['hap1'][pcol.pos+1]=(len(pileup_dict['hap1'][pcol.pos+1]),hap1_freq)
                
    
    pileup_list={'hap0':[],'hap1':[]}
    
    ref_df[['ref']]=ref_df[['ref']].applymap(lambda x:mapping[x])

    for tp in ['hap0','hap1']:
        pos_list=output[tp].keys()

        for v_pos in pos_list:

            rlist=np.array([ref_df.loc[x].ref for x in list(range(v_pos-window,v_pos+window+1)) if x in pileup_dict[tp].keys()])

            p_df=pd.DataFrame({x:pileup_dict[tp][x] for x in list(range(v_pos-window,v_pos+window+1)) if x in pileup_dict[tp].keys()})
            f_df=pd.DataFrame({x:features[tp][x] for x in list(range(v_pos-window,v_pos+window+1)) if x in features[tp].keys()})

            p_df.fillna('N',inplace=True)
            f_df.fillna(0,inplace=True)

            p_df=p_df.applymap(lambda x: mapping[x])

            if p_df.shape[0]>du_lim:
                    p_df=p_df.sample(n=du_lim,replace=False, random_state=1)

            tmp=np.array(p_df.columns())
            v_pos_range=(sum(tmp<v_pos),sum(tmp>v_pos))
            p_df=p_df.assign(f = p_df.index.map(lambda x:flag_dict[x]))
            p_df.sort_values([v_pos,'f'],inplace=True,ascending=False)
            p_df_flag=p_df.f
            p_df.drop('f', axis=1,inplace=True)
            p_df.dropna(subset=[v_pos],inplace=True)
            f_df=f_df.reindex(p_df.index)
            f_mat=np.array(f_df)
            p_mat=np.array(p_df)
            ref=(p_mat==rlist)
            p_mat=np.eye(6)[p_mat]
            p_mat=p_mat[:,:,:4]
            strand=np.array(p_df_flag)[:,np.newaxis]+np.zeros(p_mat.shape[:-1])
            data=np.dstack([p_mat,ref,f_mat,strand])
            
            cols1=p_df.columns.to_list()
            
            if data.shape[0]<du_lim:
                    tmp=np.zeros((du_lim-data.shape[0],data.shape[1],data.shape[2]))
                    data=np.vstack((data,tmp))
            
            data=np.hstack([np.zeros([32,window-sum(p_df.columns<v_pos),7]),data,np.zeros([32,window-sum(p_df.columns>v_pos),7])])
            if data.shape[1]!=17:
                print(v_pos)
                print('HAYE MAIN MAR GAYI')
                return 'asdfg'
                      
            pileup_list[tp].append((v_pos,ref_df.loc[v_pos].ref, data.astype(np.int16),output[tp][v_pos][0],output[tp][v_pos][1]))


    return pileup_list
    
    
def generate(params,mode='training'):
    cores=params['cpu']
    mode=params['mode']
    chrom=params['chrom']
    start,end=params['start'],params['end']
    
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}
    
    if mode in ['training','train']: #'pos,ref,seq,names'
        print('starting training pileups',flush=True)
        bcf_in = VariantFile(params['vcf_path'])
        tr_pos={}
        total_pos=[]
        for rec in bcf_in.fetch(chrom,start,end+1):
            total_pos.append(rec.pos)
            if len(rec.alleles[1])==1 and len(rec.ref)==1:
                if (rec.samples.items()[0][1].get('GT') in [(1,0),(0,1)] and rec.samples.items()[0][1].get('PS')) or rec.samples.items()[0][1].get('GT')==(1,1):
                    tr_pos[rec.pos]=mapping[rec.alleles[rec.samples.items()[0][1].get('GT')[0]]], mapping[rec.alleles[rec.samples.items()[0][1].get('GT')[1]]]

        tr_df=pd.DataFrame(tr_pos).transpose()
        tr_df['pos']=tr_df.index
        
        params['tr_df']=(tr_df,set(total_pos))
        
        pool = mp.Pool(processes=cores)
        fname='%s.pileups' %params['chrom']
        file_list={}
        for i in ['pos',0,15]:
            suffix='pos' if i == 'pos' else 'neg.%d' %i
            tmp_name=os.path.join(params['out_path'],'%s.%s' %(fname,suffix)) 
            file_list[i]=open(tmp_name , "w")
        
        for mbase in range(start,end,int(5e6)):
            print('starting pool:'+str(mbase),flush=True)
            t=time.time()                

            in_dict_list=[]

            for k in range(mbase,min(end,mbase+int(5e6)),100000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(k+100000,end)
                in_dict_list.append(d)
            results_dict = pool.map(get_training_candidates, in_dict_list)
                
            for result in results_dict:
                for i in ['pos',0,15]:
                    pileups=result[i]
                    if pileups:
                        for data in pileups:
                                pos,allele,ref,mat=data
                                mat=mat.reshape(-1)
                                s='%s%d%d%d%s' %((11-len(str(pos)))*'0',pos,allele,ref,''.join([(3-len(x))*' '+x for x in mat.astype('<U3')]))
                                file_list[i].write(s)
            results_dict=None
            elapsed=time.time()-t
            
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)                    

                
    elif mode in ['testing','test']:#'pos,depth,freq,ref,seq,names'

        print('starting pileups',flush=True)
        pool = mp.Pool(processes=cores)        
        out_file={}
        out_file['hap0']=open(os.path.join(params['out_path'],'%s.pileups.test.hap0' %params['chrom']) , "w")
        out_file['hap1']=open(os.path.join(params['out_path'],'%s.pileups.test.hap1' %params['chrom']) , "w")
        start,end=params['start'],params['end']
        
        for mbase in range(start,end,int(5e6)):
            print('starting pool:'+str(mbase),flush=True)
            t=time.time()                

            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(5e6)),100000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(end,k+100000)
                in_dict_list.append(d)
            results_dict = pool.map(get_testing_candidates, in_dict_list)
            
            for result in results_dict:
                for tp in ['hap0','hap1']:
                    pileups=result[tp]
                    if pileups:
                        for data in pileups:
                            pos,ref,mat,dp,freq=data
                            mat=mat.reshape(-1)
                            s='%s%d%d%s%s%d%.3f' %((11-len(str(pos)))*'0',pos,ref,''.join([(3-len(x))*' '+x for x in mat.astype('<U3')]),(3-len(str(pos)))*'0',dp,freq)
                            out_file[tp].write(s)

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
    parser.add_argument("-bed", "--bed", help="BED file")
    parser.add_argument("-hap", "--hap_path", help="Phase Info file")
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

    in_dict={'mode':args.mode, 'chrom':chrom,'start':start,'end':end,\
         'sam_path':args.bam, 'fasta_path':args.ref, 'vcf_path':args.vcf,\
             'out_path':args.output, 'window':args.window, 'depth':args.depth,\
             'cpu':args.cpu, 'bed':args.bed,'hap_path':args.hap_path,\
            'mincov':args.mincov, 'nbr_size':args.nbr_size}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
    

    
    
