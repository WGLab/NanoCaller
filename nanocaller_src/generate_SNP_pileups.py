import pysam, random
from collections import Counter
import numpy as np
from intervaltree import Interval, IntervalTree

def get_cnd_pos(v_pos,cnd_pos, seq='ont'):
    if seq=='ont':
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

    elif seq=='ul_ont':
        ls=cnd_pos[abs(cnd_pos-v_pos)<100000] 

        ls1_0= [p for p in ls if (p>=v_pos-2000) &  (p<v_pos)][:2]
        ls1_1= [p for p in ls if (p>=v_pos-5000) &  (p<v_pos-2000)][-2:]
        ls1_2= [p for p in ls if (p>=v_pos-10000) & (p<v_pos-5000)][-3:]
        ls1_3= [p for p in ls if (p>=v_pos-20000) & (p<v_pos-10000)][-3:]
        ls1_4= [p for p in ls if (p>=v_pos-40000) & (p<v_pos-20000)][-4:]
        ls1_5= [p for p in ls if (p>=v_pos-50000) & (p<v_pos-40000)][-3:]
        ls1_6= [p for p in ls if                  (p<v_pos-50000)][-3:]

        ls2_0= [p for p in ls if (p>v_pos) & (p<=v_pos+2000)][-2:]
        ls2_1= [p for p in ls if (p>v_pos+2000) & (p<=v_pos+5000)][:2]
        ls2_2= [p for p in ls if (p>v_pos+5000) & (p<=v_pos+10000)][:3]
        ls2_3= [p for p in ls if (p>v_pos+10000) & (p<=v_pos+20000)][:3]
        ls2_4= [p for p in ls if (p>v_pos+20000) & (p<=v_pos+40000)][:4]
        ls2_5= [p for p in ls if (p>v_pos+40000) & (p<=v_pos+50000)][:3]
        ls2_6= [p for p in ls if (p>v_pos+50000)][:3]

        ls_total_1=sorted(ls1_0+ls1_1+ls1_2+ls1_3+ls1_4+ls1_5+ls1_6)
        ls_total_2=sorted(ls2_0+ls2_1+ls2_2+ls2_3+ls2_4+ls2_5+ls2_6)
    
    elif seq=='ul_ont_extreme':
        ls=cnd_pos[abs(cnd_pos-v_pos)<300000] 

        ls1_0= [p for p in ls if (p>=v_pos-10000) &  (p<v_pos)][:2]
        ls1_1= [p for p in ls if (p>=v_pos-20000) & (p<v_pos-10000)][-2:]
        ls1_2= [p for p in ls if (p>=v_pos-50000) & (p<v_pos-20000)][-3:]
        ls1_3= [p for p in ls if (p>=v_pos-75000) & (p<v_pos-50000)][-3:]
        ls1_4= [p for p in ls if (p>=v_pos-100000) & (p<v_pos-75000)][-4:]
        ls1_5= [p for p in ls if (p>=v_pos-200000) & (p<v_pos-100000)][-4:]
        ls1_6= [p for p in ls if                  (p<v_pos-200000)][-2:]

        ls2_0= [p for p in ls if (p>v_pos) & (p<=v_pos+10000)][-2:]
        ls2_1= [p for p in ls if (p>v_pos+10000) & (p<=v_pos+20000)][:2]
        ls2_2= [p for p in ls if (p>v_pos+20000) & (p<=v_pos+50000)][:3]
        ls2_3= [p for p in ls if (p>v_pos+50000) & (p<=v_pos+75000)][:3]
        ls2_4= [p for p in ls if (p>v_pos+75000) & (p<=v_pos+100000)][:4]
        ls2_5= [p for p in ls if (p>v_pos+100000) & (p<=v_pos+200000)][:4]
        ls2_6= [p for p in ls if (p>v_pos+200000)][:2]

        ls_total_1=sorted(ls1_0+ls1_1+ls1_2+ls1_3+ls1_4+ls1_5+ls1_6)
        ls_total_2=sorted(ls2_0+ls2_1+ls2_2+ls2_3+ls2_4+ls2_5+ls2_6)
        
    elif seq=='pacbio':
        ls=cnd_pos[abs(cnd_pos-v_pos)<20000]
   
        ls1_0= [p for p in ls if (p>=v_pos-2000) &  (p<v_pos)][:4]
        ls1_1= [p for p in ls if (p>=v_pos-5000) &  (p<v_pos-2000)][-5:]
        ls1_2= [p for p in ls if (p>=v_pos-10000) & (p<v_pos-5000)][-5:]
        ls1_3= [p for p in ls if (p>=v_pos-20000) & (p<v_pos-10000)][-6:]
        
        ls2_0= [p for p in ls if (p>v_pos) & (p<=v_pos+2000)][-4:]
        ls2_1= [p for p in ls if (p>v_pos+2000) & (p<=v_pos+5000)][:5]
        ls2_2= [p for p in ls if (p>v_pos+5000) & (p<=v_pos+10000)][:5]
        ls2_3= [p for p in ls if (p>v_pos+10000) & (p<=v_pos+20000)][:6]
        
        ls_total_1=sorted(ls1_0+ls1_1+ls1_2+ls1_3)
        ls_total_2=sorted(ls2_0+ls2_1+ls2_2+ls2_3)

    return ls_total_1, ls_total_2

def get_snp_testing_candidates(dct, region):
    base_to_num_map={'*':4,'A':0,'G':1,'T':2,'C':3,'N':4}
    
    chrom=region['chrom']
    start=region['start']
    end=region['end']
    
    exclude_intervals=None
        
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
    
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-50000),end+50000+1),fastafile.fetch(chrom,max(1,start-50000)-1,end+50000)) }
    
    pileup_dict={}
    
    nbr_sites=[]
    
    output={}
    
    if dct['supplementary']:
        flag=0x4|0x100|0x200|0x400
    else:
        flag=0x4|0x100|0x200|0x400|0x800
        
    for pcol in samfile.pileup(chrom,max(0,start-1-50000),end+50000,min_base_quality=0,\
                                           flag_filter=flag,truncate=True):
            v_pos=pcol.pos+1
            r=ref_dict[v_pos]
            
            if not ex_bed(exclude_intervals, v_pos) and r in 'AGTC':
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                name=pcol.get_query_names()
                n=pcol.get_num_aligned()
                
                alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=r and x[0] in 'AGTC')]+[0])/n
                
                
                
                if n>=dct['mincov']:
                    # add to neighbor site list
                    if dct['threshold'][0]<=alt_freq and alt_freq<dct['threshold'][1]:
                        nbr_sites.append(v_pos)
                        pileup_dict[v_pos]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}
                    
                    # add to candidate sites list
                    if v_pos>=start and v_pos<=end and dct['min_allele_freq']<=alt_freq:
                        if v_pos not in pileup_dict:
                            pileup_dict[v_pos]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}
                        output[v_pos]=(n,alt_freq)

    nbr_sites=np.array(nbr_sites)
    
    pileup_list=[]
        
    pos_list=output.keys()
    current_depth=[]
    depth=0
    output_pos,output_ref,output_mat,output_dp,output_freq=[],[],[],[],[]
    if pos_list:

        for v_pos in pos_list:
            np.random.seed(812)
            ls_total_1, ls_total_2 = get_cnd_pos(v_pos, nbr_sites, dct['seq'])
                                
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
                        tmp_mat[i][j]=pileup_dict[nb_pos][name]
                    except KeyError:
                        pass
                        
                tmp_mat[i][len(ls_total_1)]=pileup_dict[v_pos][name]
                
                for j,nb_pos in enumerate(ls_total_2):
                    
                    try:
                        tmp_mat[i][j +len(ls_total_1)+1]=pileup_dict[nb_pos][name]
                    except KeyError:
                        pass
            for nb_pos in ls_total_1+ls_total_2:
                        rlist.append((nb_pos, base_to_num_map[ref_dict[nb_pos]]))
            
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
            data=np.hstack([np.zeros([5,nbr_size-len(ls_total_1),5]),data,np.zeros([5, nbr_size-len(ls_total_2), 5])]).astype(np.int32)
            
            output_pos.append(v_pos)
            output_ref.append(base_to_num_map[ref_dict[v_pos]])
            output_mat.append(data)
            output_dp.append(output[v_pos][0])
            output_freq.append(output[v_pos][1])
            current_depth.append(len(sample))
            
        if len(output_pos)>0:
            output_mat=np.array(output_mat).astype(np.float32)    
            output_pos=np.array(output_pos)

            output_ref=np.eye(max(4,np.max(output_ref)+1))[np.array(output_ref)].astype(np.int32)
            output_ref=output_ref[:,:4]

            output_dp=np.array(output_dp)
            output_freq=np.array(output_freq)
            depth=np.mean(current_depth)
    
    return (output_pos,output_ref,output_mat,output_dp,output_freq,depth)