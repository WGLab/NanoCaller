import sys, pysam, time, os, copy, argparse, random
from collections import Counter
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
    output_ref, output_seq={},{}
    for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            n=pcol.get_num_aligned()
            r=rlist[pcol.pos+1-start]

            if r in 'AGTC' and n>=dct['mincov']:
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=r and x[0] in 'AGTC')]+[0])/n

                if dct['threshold'][0]<=alt_freq and alt_freq<dct['threshold'][1]:
                    name=pcol.get_query_names()
                    output_ref[pcol.pos+1]=base_to_num_map[r]
                    output_seq[pcol.pos+1]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}
    return output_ref, output_seq

def get_snp_testing_candidates(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    
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
        
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    threshold=dct['threshold']
    
    nbr_size=20
    
    cnd_ref=dct['cnd_ref']
    cnd_seq=dct['cnd_seq']
    
    cnd_pos=np.array(list(cnd_ref.keys()))
    
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    pileup_dict={}
    
    output={}
    
    if dct['supplementary']:
        flag=0x4|0x100|0x200|0x400
    else:
        flag=0x4|0x100|0x200|0x400|0x800
        
    for pcol in samfile.pileup(chrom,max(0,start-1),end,min_base_quality=0,\
                                           flag_filter=flag,truncate=True):
            
            r=ref_dict[pcol.pos+1]
            if in_bed(include_intervals, pcol.pos+1) and not ex_bed(exclude_intervals, pcol.pos+1) and r in 'AGTC':
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()
                n=pcol.get_num_aligned()

                if n>=dct['mincov'] and pcol.pos+1>=start and pcol.pos+1<=end:
                    alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=r and x[0] in 'AGTC')]+[0])/n

                    if dct['min_allele_freq']<=alt_freq:
                        pileup_dict[pcol.pos+1]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}
                        output[pcol.pos+1]=(n,alt_freq)
                    
    pileup_list=[]
        
    pos_list=output.keys()

    output_pos,output_ref,output_mat,output_dp,output_freq=[],[],[],[],[]
    if pos_list:

        for v_pos in pos_list:
            np.random.seed(812)
            if dct['seq']=='ont':
                ls=cnd_pos[abs(cnd_pos-v_pos)<50000] 
                
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
            
                ls_total_1=sorted(ls1_0+ls1_1+ls1_2+ls1_3+ls1_4)
                ls_total_2=sorted(ls2_0+ls2_1+ls2_2+ls2_3+ls2_4)

            else:
                ls=cnd_pos[abs(cnd_pos-v_pos)<20000] 
                
                ls_total_1= [p for p in ls if (p>=v_pos-20000) &  (p<v_pos)][-20:]
                ls_total_2= [p for p in ls if (p>v_pos) & (p<=v_pos+20000)][:20]
                                
            nbr_dict={}

            rlist=[(v_pos,base_to_num_map[ref_dict[v_pos]]) ]
            
            sample=list(pileup_dict[v_pos].keys())
            
            if len(sample) > dct['maxcov']:
                sample=random.sample(sample,min(len(sample), dct['maxcov']))

            sample=set(sample)
            nbr_dict={}
            
            tmp_mat=np.array([[4]*(len(ls_total_1)+1+len(ls_total_2))]*len(sample))
            
            for i,name in enumerate(sample):
                for j,nb_pos in enumerate(ls_total_1):
                    try:    
                        tmp_mat[i][j]=cnd_seq[nb_pos][name]
                    except KeyError:
                        pass
                        
                tmp_mat[i][len(ls_total_1)]=pileup_dict[v_pos][name]
                
                for j,nb_pos in enumerate(ls_total_2):
                    
                    try:
                        tmp_mat[i][j +len(ls_total_1)+1]=cnd_seq[nb_pos][name]
                    except KeyError:
                        pass
            for nb_pos in ls_total_1+ls_total_2:
                        rlist.append((nb_pos, cnd_ref[nb_pos]))
            
            rlist=sorted(rlist, key=lambda x:x[0])
            total_rlist=np.array([x[1] for x in rlist])
            
            if len(total_rlist)<dct['min_nbr_sites']:
                continue
            
            mat=np.dstack([np.sum(np.eye(5)[tmp_mat[tmp_mat[:,len(ls_total_1)]==i]],axis=0) for i in range(4)]).transpose(2,0,1)[:,:,:4]
            
            total_ref=np.eye(5)[total_rlist.astype(int)]
            total_ref[:,4]=0
            total_ref=total_ref[np.newaxis,:]
            mat=np.dstack([mat,np.zeros([4,mat.shape[1]])+np.eye(4)[base_to_num_map[ref_dict[v_pos]]][:,np.newaxis]])
            data=np.vstack([total_ref,np.multiply(mat,1-2*total_ref)])
            data=np.hstack([np.zeros([5,nbr_size-len(ls_total_1),5]),data,np.zeros([5, nbr_size-len(ls_total_2), 5])]).astype(np.int8)
            
            output_pos.append(v_pos)
            output_ref.append(base_to_num_map[ref_dict[v_pos]])
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

def generate(params, pool):
    chrom=params['chrom']
    threshold=params['threshold']
    start,end=params['start'],params['end']    
    
    print('generating SNP pileups for %s:%d-%d' %(chrom,start,end) ,flush=True)
    
    if params['include_bed']:
        tbx = pysam.TabixFile(params['include_bed'])
        bed_intervals=IntervalTree(Interval(int(row[1]), int(row[2]), "%s" % (row[1])) for row in tbx.fetch(chrom, parser=pysam.asBed()))
        bed_intervals=IntervalTree(bed_intervals.overlap(start,end))
        
        if not bed_intervals:
            return [],[],[],[],[]
    
        else:
            start=max(start, min(x[0] for x in bed_intervals))
            end=min(end, max(x[1] for x in bed_intervals))
            
            
    
    in_dict_list=[]
    for k in range(max(1,start-50000),end+50000,100000):
        d = copy.deepcopy(params)
        d['start']=k
        d['end']=min(end+50000,k+100000)
        in_dict_list.append(d)
        
    nbr_results = pool.map(get_nbr, in_dict_list)
    
    total_nbr_ref, total_nbr_seq={}, {}
    for pair in nbr_results:
        total_nbr_ref.update(pair[0])
        total_nbr_seq.update(pair[1])
    
    params['cnd_ref']=total_nbr_ref
    params['cnd_seq']=total_nbr_seq
        
    in_dict_list=[]
    for k in range(start,end,100000):
        d = copy.deepcopy(params)
        d['start']=k
        d['end']=min(end,k+100000)
        in_dict_list.append(d)
    results = pool.map(get_snp_testing_candidates, in_dict_list)

    if sum([len(res[0]) for res in results])==0:
        return [],[],[],[],[]
    
    pos=np.vstack([res[0][:,np.newaxis] for res in results if len(res[0])>0])
    mat=np.vstack([res[2] for res in results if len(res[0])>0])
    
    ref=np.vstack([res[1] for res in results if len(res[0])>0 ])
    dp=np.vstack([res[3][:,np.newaxis] for res in results if len(res[0])>0])

    freq=np.vstack([res[4][:,np.newaxis] for res in results if len(res[0])>0])
    
    return pos,mat,ref,dp,freq