import sys,pysam, time,os,copy,argparse,subprocess, random, re, datetime
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from subprocess import Popen, PIPE, STDOUT
from intervaltree import Interval, IntervalTree
import collections
import parasail

sub_mat=parasail.matrix_create('AGTC',20,-10)

def msa(seq_list, ref, v_pos, mincov, maxcov):
    mapping={'A':0,'G':1,'T':2,'C':3,'-':4}
    rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}
    
    np.random.seed(812)
    sample=list(seq_list.keys())
            
    if len(sample) > maxcov:
        sample=random.sample(sample,min(len(sample),maxcov))

    sample=sorted(sample)

    fa_tmp_file=''.join(['>%s_SEQ\n%s\n'%(read_name,seq_list[read_name]) for read_name in sample])


    fa_tmp_file+='>ref_SEQ\n%s' %ref
    
    gap_penalty=1.0
    msa_process =Popen(['muscle', '-quiet','-gapopen','%.1f' %gap_penalty,'-maxiters', '1' ,'-diags1'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    hap_file=msa_process.communicate(input=fa_tmp_file.encode('utf-8'))

    if len(hap_file)==0:
        return (0,0,None,None,None)


    tmp=hap_file[0].decode('utf-8')[1:].replace('\n','').split('>')

    zz_0=[]
    for seq in tmp:
        p1,p2=seq.split('_SEQ')
        if p1!='ref':
            zz_0.append(p2)
        else:
            ref_real_0=p2

    if len(zz_0)<mincov:
        return (0,0,None,None,None)

    cnt_del=0
    
    try:
        ref_real_0_mat=np.eye(5)[[mapping[x] for x in ref_real_0]]
    except UnboundLocalError:
        return (0,0,None,None,None)
    
    mat=np.array([[mapping[c] for c in x] for x in zz_0])
    h0_mat=np.sum(np.eye(5)[mat],axis=0).astype(np.float32)
    alt_mat=(h0_mat/(np.sum(h0_mat,axis=1)[:,np.newaxis]))
    
    tmp_mat=np.copy(alt_mat)
    tmp_mat[:,4]=tmp_mat[:,4]-0.01
    
    cns=''.join([rev_base_map[x] for x in np.argmax(tmp_mat,axis=1)]).replace('-','')
    ref_seq=ref_real_0.replace('-','')
    
    alt_mat-=ref_real_0_mat
    
    final_mat=np.dstack([alt_mat, ref_real_0_mat])[:128,:,:].transpose(1,0,2)
    if final_mat.shape[1]<128:
        final_mat=np.hstack((final_mat,np.zeros((5, 128-final_mat.shape[1], 2))))
        
    return (1,1,final_mat,cns,ref_seq)



def allele_prediction(alt, ref_seq, max_range):
    global sub_mat
    cigar=parasail.nw_trace(alt, ref_seq, 9, 1, sub_mat).cigar
    cigar_op=[(x & 0xf, x >> 4) for x in cigar.seq]    
    
    indel=False
    ref_cnt=[0]*10
    alt_cnt=[0]*10
    
    mis_match_cnt_before_indel=False
    mis_match_cnt_after_indel=(0, 0)
    for op, cnt in cigar_op:
        if op==8 or op==7:
            ref_cnt[op]+=cnt
            alt_cnt[op]+=cnt
            
            if indel:
                mis_match_cnt_after_indel[op-7]+=cnt
            else:
                mis_match_cnt_before_indel=True
            
        if op==1:
            alt_cnt[op]+=cnt
            mis_match_cnt_after_indel=[0, 0]
            indel=True
        
        if op==2:
            ref_cnt[op]+=cnt
            mis_match_cnt_after_indel=[0, 0]
            indel=True
        
        if indel==False and sum(ref_cnt)>=max_range+10:
            if ref_cnt[8]:
                out_len=sum(ref_cnt) if op==8 else sum(ref_cnt)-cnt
                return ref_seq[:out_len], alt[:out_len]
            else:
                return (None, None)
            
            
        if indel==True:
            if sum(mis_match_cnt_after_indel)>20:
                break
                
    ref_out_len=sum(ref_cnt) if op==8 else sum(ref_cnt)-cnt
    alt_out_len=sum(alt_cnt) if op==8 else sum(alt_cnt)-cnt
    
    if not mis_match_cnt_before_indel:
        ref_out_len+=1
        alt_out_len+=1
        
    return ref_seq[:ref_out_len], alt[:alt_out_len]

def get_indel_testing_candidates(dct, chunk):
    mapping={'A':0,'G':1,'T':2,'C':3,'-':4}
    rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

    init_time=time.time()
    start_time=str(datetime.datetime.now())
    
    window_before,window_after=0,160
    
    if dct['seq']=='pacbio':
        window_after=260
    
    chrom=chunk['chrom']
    start=chunk['start']
    end=chunk['end']
    
    sam_path=chunk['sam_path']
    fasta_path=dct['fasta_path']
    samfile = pysam.Samfile(sam_path, reference_filename=fasta_path)
    fastafile=pysam.FastaFile(fasta_path)

    window_size=dct['win_size']
    small_size=dct['small_win_size']
        
    mincov,maxcov=dct['mincov'],dct['maxcov']
    ins_t,del_t=dct['ins_t'],dct['del_t']
    
    exclude_intervals = None
    
    
    if dct['exclude_bed']:
        tbx = pysam.TabixFile(dct['exclude_bed'])
        try:
            exclude_intervals=IntervalTree(Interval(int(row[1]), int(row[2]), "%s" % (row[1])) for row in tbx.fetch(chrom, parser=pysam.asBed()))
        
            def ex_bed(tree, pos):
                return tree.overlaps(pos)
        
        except ValueError:
            def ex_bed(tree, pos):
                return False 
    else:
        def ex_bed(tree, pos):
            return False 
        
    ref_dict={j:s.upper() if s in 'AGTC' else 'N' for j,s in zip(range(max(1,start-200),end+400+1),fastafile.fetch(chrom,max(1,start-200)-1,end+400)) }
    
    chrom_length=fastafile.get_reference_length(chrom)
    
    hap_dict={1:[],2:[]}
    phase_dict={}
    for pread in samfile.fetch(chrom, max(0, start-100000), end + 1000, multiple_iterators=False):
        if pread.has_tag('HP'):
            phase_dict[pread.qname]=pread.get_tag('PS')
            hap_dict[pread.get_tag('HP')].append(pread.qname)
        else:
            phase_dict[pread.qname]=None

    hap_reads_0=set(hap_dict[1])
    hap_reads_1=set(hap_dict[2])
    
    if dct['supplementary']:
        flag=0x4|0x100|0x200|0x400
    else:
        flag=0x4|0x100|0x200|0x400|0x800
        
    output_pos,output_data_0,output_data_1,output_data_total,alleles, output_phase=[], [], [], [], [], []
    
    del_queue_0,del_queue_1=collections.deque(window_size*[set()], window_size),collections.deque(window_size*[set()], window_size)
    ins_queue_0,ins_queue_1=collections.deque(window_size*[set()], window_size),collections.deque(window_size*[set()], window_size)
    position_queue=collections.deque(window_size*[{}], window_size)
    
    
    del_queue_small_0,del_queue_small_1=collections.deque(small_size*[set()], small_size),collections.deque(small_size*[set()], small_size)
    ins_queue_small_0,ins_queue_small_1=collections.deque(small_size*[set()], small_size),collections.deque(small_size*[set()], small_size)
    position_queue_small=collections.deque(small_size*[{}], small_size)

    variants={}
    extra_variants={}
    
    max_range={0:max(10,window_size),1:10}
    
    count=0
    prev=0
    for pcol in samfile.pileup(chrom,max(0,start-1),end,min_base_quality=0,\
                                           flag_filter=flag,truncate=True, multiple_iterators=False):
            v_pos=pcol.pos+1
                
            if not ex_bed(exclude_intervals, v_pos):                
                read_names=pcol.get_query_names()
                read_names_0=set(read_names) & hap_reads_0
                read_names_1=set(read_names) & hap_reads_1
                len_seq_tot=len(read_names)
                len_seq_0=len(read_names_0)
                len_seq_1=len(read_names_1)

                ins_set_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))>2  and int(''.join(filter(str.isdigit, s)))<=50}
                ins_set_small_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))<=10}

                del_set_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))>2 and int(''.join(filter(str.isdigit, s)))<=50}
                del_set_small_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))<=10}

                ins_set_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))>2 and int(''.join(filter(str.isdigit, s)))<=50}
                ins_set_small_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))<=10}

                del_set_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))>2 and int(''.join(filter(str.isdigit, s)))<=50}
                del_set_small_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))<=10}

                del_queue_0.append(del_set_0)
                del_queue_1.append(del_set_1)

                ins_queue_0.append(ins_set_0)
                ins_queue_1.append(ins_set_1)

                del_queue_small_0.append(del_set_small_0)
                del_queue_small_1.append(del_set_small_1)

                ins_queue_small_0.append(ins_set_small_0)
                ins_queue_small_1.append(ins_set_small_1)

                if v_pos<=prev:
                    continue
                    
                if len_seq_0>=mincov  and len_seq_1>=mincov:

                    del_freq_0=len(set.union(*del_queue_0))/len_seq_0 if len_seq_0>0 else 0
                    ins_freq_0=len(set.union(*ins_queue_0))/len_seq_0 if len_seq_0>0 else 0

                    del_freq_1=len(set.union(*del_queue_1))/len_seq_1 if len_seq_1>0 else 0
                    ins_freq_1=len(set.union(*ins_queue_1))/len_seq_1 if len_seq_1>0 else 0

                    del_freq_small_0=len(set.union(*del_queue_small_0))/len_seq_0 if len_seq_0>0 else 0
                    ins_freq_small_0=len(set.union(*ins_queue_small_0))/len_seq_0 if len_seq_0>0 else 0

                    del_freq_small_1=len(set.union(*del_queue_small_1))/len_seq_1 if len_seq_1>0 else 0
                    ins_freq_small_1=len(set.union(*ins_queue_small_1))/len_seq_1 if len_seq_1>0 else 0

                    if max([del_freq_0,del_freq_1])>=del_t or max([ins_freq_0,ins_freq_1])>=ins_t:
                        prev=v_pos+window_size
                        variants[max(1,v_pos-window_size)]=0
                        count+=1

                    elif max([del_freq_small_0,del_freq_small_1])>=del_t or max([ins_freq_small_0,ins_freq_small_1])>=ins_t or (del_freq_small_0+ins_freq_small_0)>=0.9 or (del_freq_small_1+ins_freq_small_1)>=0.9:

                        prev=v_pos+10
                        variants[max(1,v_pos-10)]=1
                        count+=1
                
            
                elif dct['impute_indel_phase'] and len_seq_tot >=2*mincov:
                        seq_v2=[x.upper() for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)]
                        seq=[x[:2] for x in seq_v2]
                        seq_tot=''.join(seq)

                        del_freq_tot=(seq_tot.count('-')+seq_tot.count('*'))/len_seq_tot if len_seq_tot>0 else 0
                        ins_freq_tot=seq_tot.count('+')/len_seq_tot if len_seq_tot>0 else 0
                        if (del_t<=del_freq_tot or ins_t<=ins_freq_tot):
                            groups={}
                            for s,n in zip(seq_v2,read_names):
                                if s not in groups:
                                    groups[s]=[]
                                groups[s].append(n)

                            counts=sorted([(x,len(groups[x])) for x in groups],key=lambda x:x[1],reverse=True)
                            
                            if counts[0][1]<=0.8*len_seq_tot:
                                read_names_0=set(groups[counts[0][0]])
                                read_names_1=set(groups[counts[1][0]]) if counts[1][1]>=mincov else set(read_names)-read_names_0
                            else:
                                read_names_0=groups[counts[0][0]][:counts[0][1]//2]
                                read_names_1=groups[counts[0][0]][counts[0][1]//2:]
                            if len(read_names_0)>=mincov and len(read_names_1)>=mincov:
                                prev=v_pos+10
                                variants[max(1,v_pos-10)]=1
                                extra_variants[max(1,v_pos-10)]=(read_names_0,read_names_1)
                                count+=1

    for pcol in samfile.pileup(chrom,max(0,start-10-window_size),end,min_base_quality=0, flag_filter=flag,truncate=True, multiple_iterators=False):
                    
        v_pos=pcol.pos+1
        if v_pos in extra_variants:
            read_names=pcol.get_query_names()
            read_names_0, read_names_1 =extra_variants[v_pos]
            
        elif v_pos in variants:
            read_names=pcol.get_query_names()
            
            read_names_0=set(read_names) & hap_reads_0
            read_names_1=set(read_names) & hap_reads_1
            
        else:
            continue

        d={'hap0':{},'hap1':{}}
        d_tot={}

        ref=''.join([ref_dict[p] for p in range(v_pos-window_before,min(chrom_length,v_pos+window_after+1))])

        if 'N' in ref:
            continue
            
        for pread in pcol.pileups:
            dt=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]                    
            d_tot[pread.alignment.qname]=dt

            if pread.alignment.qname in read_names_0:
                d['hap0'][pread.alignment.qname]=dt

            elif pread.alignment.qname in read_names_1:
                d['hap1'][pread.alignment.qname]=dt    


        seq_list=d['hap0']
        flag0,_, data_0,alt_0,ref_seq_0=msa(seq_list,ref,v_pos,2,dct['maxcov'])

        seq_list=d['hap1']
        flag1,_, data_1,alt_1,ref_seq_1=msa(seq_list,ref,v_pos,2,dct['maxcov'])

        seq_list = d_tot
        flag_total,indel_flag_total,data_total,alt_total,ref_seq_total=msa(seq_list,ref,v_pos,dct['mincov'],dct['maxcov'])

        if flag0 and flag1 and flag_total:
            output_pos.append(v_pos)
            output_data_0.append(data_0)
            output_data_1.append(data_1)
            output_data_total.append(data_total)
            output_phase.append(phase_dict[next(iter(d['hap0'].keys()))])

            tp=max_range[variants[v_pos]]

            alleles.append([allele_prediction(alt_0, ref_seq_0, max_range[variants[v_pos]]),\
                            allele_prediction(alt_1, ref_seq_1, max_range[variants[v_pos]]), \
                            allele_prediction(alt_total, ref_seq_total, max_range[variants[v_pos]])])

    
    if len(output_pos)==0:
        return (output_pos, output_data_0, output_data_1, output_data_total, alleles, output_phase)
    
    output_data_0=np.array(output_data_0)
    output_data_1=np.array(output_data_1)
    output_data_total=np.array(output_data_total)
        
    return (output_pos, output_data_0, output_data_1, output_data_total, alleles, output_phase)