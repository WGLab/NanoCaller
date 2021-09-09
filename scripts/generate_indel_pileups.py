import sys,pysam, time,os,copy,argparse,subprocess, random, re, datetime
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from subprocess import Popen, PIPE, STDOUT
from Bio import pairwise2
from intervaltree import Interval, IntervalTree
import collections



def pairwise(x,y):
    alignments = pairwise2.align.globalms(x, y, 2, -1.0, -0.9, -0.1)

    return alignments

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
    
    cns=''.join([rev_base_map[x] for x in np.argmax(tmp_mat,axis=1)])
    ref_seq=ref_real_0.replace('-','')
    
    alt_mat-=ref_real_0_mat
    final_mat=np.dstack([alt_mat, ref_real_0_mat])[:128,:,:]    
    return (1,1,final_mat.transpose(1,0,2),cns,ref_seq)


def allele_prediction(alt,ref_seq,max_range):
    res=pairwise(alt.replace('-',''),ref_seq.replace('-',''))
    try:
        flag=False
        best=None
        score=1e6
        for pair in res:
            current=sum([m.start() for m in re.finditer('-', pair[0])])+sum([m.start() for m in re.finditer('-', pair[1])])
            if current <score:
                best=pair
                score=current

        if best==None:
            return (None,None)
        
        ref_count=0
        alt_new=best[0]
        ref_new=best[1]
        for i in range(len(ref_new)-2):
            if ref_new[i]!='-':
                ref_count+=1
            if ref_count>=max_range+10 and flag==False:
                break
            if ref_new[i]=='-' or alt_new[i]=='-':
                flag=True
            if ref_new[i:i+10]==alt_new[i:i+10] and flag and ref_count>=max_range+10:
                break

        i+=1
        s2=alt_new[:i].replace('-','')
        s1=ref_new[:i].replace('-','')
    
    
        s1=s1[0]+'.'+s1[1:]+'|'
        s2=s2[0]+'.'+s2[1:]+'|'
    except IndexError:
        return (None,None)
    if s1==s2:
        return (None,None)
    i=-1
        
    while s1[i]==s2[i] and s1[i]!='.' and s2[i]!='.':
        i-=1

    ref_out=s1[:i+1].replace('.','')
    allele_out=s2[:i+1].replace('.','')
    
    return (ref_out,allele_out)

def get_indel_testing_candidates(dct):
    mapping={'A':0,'G':1,'T':2,'C':3,'-':4}
    rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

    init_time=time.time()
    start_time=str(datetime.datetime.now())
    
    window_before,window_after=0,160
    
    if dct['seq']=='pacbio':
        window_after=260
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    window_size=dct['win_size']
    small_size=dct['small_win_size']
        
    mincov,maxcov=dct['mincov'],dct['maxcov']
    ins_t,del_t=dct['ins_t'],dct['del_t']
    
    include_intervals, exclude_intervals=None, None
    
    if dct['include_bed']:
        tbx = pysam.TabixFile(dct['include_bed'])
        include_intervals=IntervalTree(Interval(int(row[1]), int(row[2]), "%s" % (row[1])) for row in tbx.fetch(chrom, parser=pysam.asBed()))
        
        def in_bed(tree,pos):            
            return tree.overlaps(pos)
        
        include_intervals=IntervalTree(include_intervals.overlap(start,end))
        
        if not include_intervals:
            return [],[],[],[],[]
    
        else:
            start=max(start, min(x[0] for x in include_intervals))
            end=min(end, max(x[1] for x in include_intervals))
        
    else:
        def in_bed(tree, pos):
            return True
    
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
        
    ref_dict={j:s.upper() if s in 'AGTC' else '' for j,s in zip(range(max(1,start-200),end+400+1),fastafile.fetch(chrom,max(1,start-200)-1,end+400)) }
    
    hap_dict={1:[],2:[]}
    
    for pread in samfile.fetch(chrom, max(0,dct['start']-100), dct['end']+100):
        if pread.has_tag('HP'):
            hap_dict[pread.get_tag('HP')].append(pread.qname)

    hap_reads_0=set(hap_dict[1])
    hap_reads_1=set(hap_dict[2])
    
    if dct['supplementary']:
        flag=0x4|0x100|0x200|0x400
    else:
        flag=0x4|0x100|0x200|0x400|0x800
        
    output_pos,output_data_0,output_data_1,output_data_total,alleles=[],[],[],[],[]
    
    del_queue_0,del_queue_1=collections.deque(window_size*[set()], window_size),collections.deque(window_size*[set()], window_size)
    ins_queue_0,ins_queue_1=collections.deque(window_size*[set()], window_size),collections.deque(window_size*[set()], window_size)
    position_queue=collections.deque(window_size*[{}], window_size)
    
    
    del_queue_small_0,del_queue_small_1=collections.deque(small_size*[set()], small_size),collections.deque(small_size*[set()], small_size)
    ins_queue_small_0,ins_queue_small_1=collections.deque(small_size*[set()], small_size),collections.deque(small_size*[set()], small_size)
    position_queue_small=collections.deque(small_size*[{}], small_size)

    variants={}
    max_range={0:max(10,window_size),1:10}
    
    count=0
    prev=0
    for pcol in samfile.pileup(chrom,max(0,start-1),end,min_base_quality=0,\
                                           flag_filter=flag,truncate=True):
            v_pos=pcol.pos+1
            if count>200:
                print('%s: Skipping region %s:%d-%d due to poor alignment.' %(str(datetime.datetime.now()), chrom, v_pos, end), flush =True)
                break
                
            if in_bed(include_intervals, v_pos) and not ex_bed(exclude_intervals, v_pos):                
                read_names=pcol.get_query_names()
                read_names_0=set(read_names) & hap_reads_0
                read_names_1=set(read_names) & hap_reads_1
                len_seq_tot=len(read_names)
                len_seq_0=len(read_names_0)
                len_seq_1=len(read_names_1)

                ins_set_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))>2}
                ins_set_small_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))<=10}

                del_set_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))>2}
                del_set_small_0={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_0 and int(''.join(filter(str.isdigit, s)))<=10}

                ins_set_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))>2}
                ins_set_small_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '+' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))<=10}

                del_set_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))>2}
                del_set_small_1={n for n,s in zip(read_names,pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)) if '-' in s and n in read_names_1 and int(''.join(filter(str.isdigit, s)))<=10}

                del_queue_0.append(del_set_0)
                del_queue_1.append(del_set_1)

                ins_queue_0.append(ins_set_0)
                ins_queue_1.append(ins_set_1)

                del_queue_small_0.append(del_set_small_0)
                del_queue_small_1.append(del_set_small_1)

                ins_queue_small_0.append(ins_set_small_0)
                ins_queue_small_1.append(ins_set_small_1)


                if v_pos>prev and len_seq_0>=mincov  and len_seq_1>=mincov:

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
            
                    elif dct['seq']=='pacbio' and len_seq_tot >=2*mincov:
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
                                read_names_0=groups[counts[0][0]]
                                read_names_1=groups[counts[1][0]]
                            else:
                                read_names_0=groups[counts[0][0]][:counts[0][1]//2]
                                read_names_1=groups[counts[0][0]][counts[0][1]//2:]

                            if len(read_names_0)>=mincov and len(read_names_1)>=mincov:
                                prev=v_pos+10
                                variants[max(1,v_pos-10)]=1
                                count+=1

    for pcol in samfile.pileup(chrom,max(0,start-10-window_size),end,min_base_quality=0, flag_filter=flag,truncate=True):
                    
        v_pos=pcol.pos+1
        
        if v_pos in variants:
            read_names=pcol.get_query_names()
            read_names_0=set(read_names) & hap_reads_0
            read_names_1=set(read_names) & hap_reads_1
           

            d={'hap0':{},'hap1':{}}
            d_tot={}
                        
            ref=''.join([ref_dict[p] for p in range(v_pos-window_before,min(end,v_pos+window_after+1))])
            
                    
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
                
                tp=max_range[variants[v_pos]]
                
                alleles.append([allele_prediction(alt_0, ref_seq_0, max_range[variants[v_pos]]),\
                                allele_prediction(alt_1, ref_seq_1, max_range[variants[v_pos]]), \
                                allele_prediction(alt_total, ref_seq_total, max_range[variants[v_pos]])])
                
                

    if len(output_pos)==0:
        return (output_pos,output_data_0,output_data_1,output_data_total,alleles)

    output_pos=np.array(output_pos)
    output_data_0=np.array(output_data_0)
    output_data_1=np.array(output_data_1)
    output_data_total=np.array(output_data_total)
    
    return (output_pos,output_data_0,output_data_1,output_data_total,alleles)