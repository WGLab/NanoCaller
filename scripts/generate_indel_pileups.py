import sys,pysam, time,os,re,copy,argparse,gzip,itertools,subprocess,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree
from tempfile import mkstemp
from subprocess import Popen, PIPE, STDOUT

os.environ["PATH"] += os.pathsep + '/home/ahsanm/lib/tcoffee/t_coffee/src'

window_before,window_after=0,160

gt_map={(0,0):0, (1,1):1, (0,1):2, (1,0):2,(1,2):3, (2,1):3}

mapping={'A':0,'G':1,'T':2,'C':3,'-':4}
revmapping={0:'A',1:'G',2:'T',3:'C',4:'-'}

allele_map={('I','I'):0, ('D','D'):1, ('N','I'):2, ('I','N'):3, ('N','D'):4, ('D','N'):5, \
           ('I','D'):6, ('D','I'):7, ('N','N'):8}

def v_type(ref,allele1,allele2):
    allele1_type,allele2_type='N','N'
    
    if len(ref)>len(allele1):
        allele1_type='D'
        
    elif len(ref)<len(allele1):
        allele1_type='I'
        
    if len(ref)>len(allele2):
        allele2_type='D'
        
    elif len(ref)<len(allele2):
        allele2_type='I'
        
    return allele_map[(allele1_type,allele2_type)]

def msa(seq_list, ref, v_pos, mincov):
    fa_tmp_file=''
    for seq_count,b in seq_list.items():
        fa_tmp_file+='>%s_SEQ\n'%seq_count
        fa_tmp_file+= '%s\n' %b


    fa_tmp_file+='>ref_SEQ\n'
    fa_tmp_file+= '%s' %ref
    
    gap_penalty=1.0
    msa_process =Popen(['muscle', '-quiet','-gapopen','%.1f' %gap_penalty,'-maxiters', '1' ,'-diags1'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    hap_file=msa_process.communicate(input=fa_tmp_file.encode('utf-8'))

    if len(hap_file)==0:
        print('hapfile length 0')


    tmp=hap_file[0].decode('utf-8')[1:].replace('\n','').split('>')

    zz_0=[]
    for seq in tmp:
        p1,p2=seq.split('_SEQ')
        if p1!='ref':
            zz_0.append(p2[:128])
        else:
            ref_real_0=p2

    if len(zz_0)<mincov:
        return (0,0,None)

    cnt_del=0
    
    try:
        ref_real_0_mat=np.eye(5)[[mapping[x] for x in ref_real_0[:128]]].transpose()
    except UnboundLocalError:
        return (0,0,None)
    
    mat=np.array([[mapping[c] for c in x] for x in zz_0])
    h0_mat=np.sum(np.eye(5)[mat],axis=0).transpose()

    return (1,1,np.dstack([h0_mat, ref_real_0_mat]))
    

def get_testing_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
            
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
    
    hap_dict={1:[],2:[]}
    
    for pread in samfile.fetch(chrom, max(0,dct['start']-100), dct['end']+100):
        if pread.has_tag('HP'):
            hap_dict[pread.get_tag('HP')].append(pread.qname)

    hap_reads_0=set(hap_dict[1])
    hap_reads_1=set(hap_dict[2])
    
    
    pileup_list=[]
    output_pos,output_data_0,output_data_1,output_data_total=[],[],[],[]
    prev=0
    for pcol in samfile.pileup(chrom,max(0,start-1),end,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            
            if pcol.pos+1>prev+1:                
                read_names=pcol.get_query_names()
                read_names_0=[x for x in read_names if x in hap_reads_0]
                read_names_1=[x for x in read_names if x in hap_reads_1]

                if len(read_names)>=dct['mincov']:
                    v_pos=pcol.pos+1
                    seq=[x[:2].upper() for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)]
                    
                    len_seq=len(seq)
                    tmp_seq=''.join(seq)
                    
                    tmp_0=[s for n,s in zip(read_names,seq) if n in read_names_0]
                    len_seq_0=len(tmp_0)
                    tmp_seq_0=''.join(tmp_0)
                    
                    tmp_1=[s for n,s in zip(read_names,seq) if n in read_names_1]
                    len_seq_1=len(tmp_1)
                    tmp_seq_1=''.join(tmp_1)
                                        
                    del_freq=(tmp_seq.count('-')+tmp_seq.count('*'))/len_seq if len_seq>0 else 0
                    ins_freq=tmp_seq.count('+')/len_seq if len_seq>0 else 0
                    

                    del_freq_0=(tmp_seq_0.count('-')+tmp_seq_0.count('*'))/len_seq_0 if len_seq_0>0 else 0
                    ins_freq_0=tmp_seq_0.count('+')/len_seq_0 if len_seq_0>0 else 0

                    del_freq_1=(tmp_seq_1.count('-')+tmp_seq_1.count('*'))/len_seq_1 if len_seq_1>0 else 0
                    ins_freq_1=tmp_seq_1.count('+')/len_seq_1 if len_seq_1>0 else 0
                    
                    ref=''.join(ref_df.reindex(range(v_pos-window_before,v_pos+window_after+1)).fillna('').ref.to_list())
                    
                    if len(ref)<128:
                        continue
                    if len_seq_0>=dct['mincov']  and len_seq_1>=dct['mincov']  and (dct['del_t']<=del_freq_0 or dct['del_t']<=del_freq_1 or dct['ins_t']<=ins_freq_0 or dct['ins_t']<=ins_freq_1):
                                            
                        
                        d={'hap0':{},'hap1':{}}
                        d_tot={}
                        for pread in pcol.pileups:
                            d_tot[pread.alignment.qname]=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]
                            if pread.alignment.qname in read_names_0:
                                d['hap0'][pread.alignment.qname]=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]

                            elif pread.alignment.qname in read_names_1:
                                d['hap1'][pread.alignment.qname]=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]
                        
                        seq_list=d['hap0']
                        flag0, indel_flag0, data_0=msa(seq_list,ref,v_pos,2)

                        
                        seq_list=d['hap1']
                        flag1,indel_flag1,data_1=msa(seq_list,ref,v_pos,2)
                        
                        seq_list = d_tot
                        flag_total,indel_flag_total,data_total=msa(seq_list,ref,v_pos,dct['mincov'])

                        if flag0 and flag1 and flag_total:
                            prev=pcol.pos+1
                            output_pos.append(v_pos)
                            output_data_0.append(data_0)
                            output_data_1.append(data_1)
                            output_data_total.append(data_total)
                            
    if len(output_pos)==0:
        return (output_pos,output_data_0,output_data_1,output_data_total)
    output_pos=np.array(output_pos)
    output_data_0=np.array(output_data_0)
    output_data_1=np.array(output_data_1)
    output_data_total=np.array(output_data_total)
    
    return (output_pos,output_data_0,output_data_1,output_data_total)
    
def generate(params,pool):
    cores=params['cpu']
    chrom=params['chrom']
    start,end=params['start'],params['end']
    
    start,end=params['start'],params['end']

    print('generating indel pileups for %s:%d-%d' %(chrom,start,end),flush=True)
    
    in_dict_list=[]
    
    for k in range(start,end,100000):
        d = copy.deepcopy(params)
        d['start']=k
        d['end']=min(end,k+100000)
        in_dict_list.append(d)
    
    results = pool.map(get_testing_candidates, in_dict_list)
    
    if sum([len(res[0]) for res in results])==0:
        return [],[],[],[]
    pos=np.vstack([res[0][:,np.newaxis] for res in results if len(res[0])>0])
    mat_0=np.vstack([res[1] for res in results if len(res[0])>0])
    mat_1=np.vstack([res[2] for res in results if len(res[0])>0])
    mat_2=np.vstack([res[3] for res in results if len(res[0])>0])
    
    return pos, mat_0, mat_1, mat_2