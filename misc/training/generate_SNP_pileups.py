import sys, pysam, time, os, copy, argparse, random, datetime
from collections import Counter
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree

base_to_num_map={'*':4,'A':0,'G':1,'T':2,'C':3,'N':4}

def in_bed(tree,pos):            
        return tree.overlaps(pos)
    
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
        
    elif seq=='new_pcb':
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
        
    elif seq=='old_pcb':
        ls=cnd_pos[abs(cnd_pos-v_pos)<20000] 

        ls_total_1= [p for p in ls if (p>=v_pos-20000) &  (p<v_pos)][-20:]
        ls_total_2= [p for p in ls if (p>v_pos) & (p<=v_pos+20000)][:20]
    
    return ls_total_1, ls_total_2
    
def get_nbr(dct, nbr_type='freq'):
    chrom=dct['chrom']
    start=max(dct['start']-50000,1)
    end=dct['end']+50000

    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    tbx = pysam.TabixFile(dct['exclude_bed'])
    exclude_intervals=IntervalTree(Interval(int(row[1]), int(row[2]), "%s" % (row[1])) for row in tbx.fetch(chrom, parser=pysam.asBed()))
    
    
    output_seq={}
    
    flag=0x4|0x100|0x200|0x400|0x800
        
    if nbr_type=='freq':
        rlist=[s for s in fastafile.fetch(chrom,start-1,end-1)]
        
        for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=flag,truncate=True):
                
                r=rlist[pcol.pos+1-start]                
                if r in 'AGTC' and not in_bed(exclude_intervals,pcol.pos+1):
                    n=pcol.get_num_aligned()
                    seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                    alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=r and x[0] in 'AGTC')]+[0])/n

                    if dct['threshold'][0]<=alt_freq and alt_freq<dct['threshold'][1] and n>=dct['mincov']:
                        name=pcol.get_query_names()
                        output_seq[pcol.pos+1]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}
                        output_seq[pcol.pos+1]['ref']=base_to_num_map[r]
    
    elif nbr_type=='gtruth':
        gt_map={(0,0):0, (1,1):0, (2,2):0, (1,2):1, (2,1):1, (0,1):1, (1,0):1, (0,2):1, (2,0):1, (1,None):0,(None,1):0}
        
        bcf_in = VariantFile(dct['vcf_path'])
        
        ground_truth={}
        for rec in bcf_in.fetch(chrom,start,end+1):
            gt=rec.samples.items()[0][1].get('GT')
            if gt_map[gt]:
                if base_to_num_map[rec.alleles[gt[0]]]<4 and base_to_num_map[rec.alleles[gt[1]]]<4:
                    ground_truth[rec.pos]=rec.ref
                         
        for pcol in samfile.pileup(chrom,start-1,end-1,min_base_quality=0, flag_filter=flag,truncate=True):
            
            if pcol.pos+1 in ground_truth:
                n=pcol.get_num_aligned()
                r=ground_truth[pcol.pos+1]


                if n>=dct['mincov']:
                    name=pcol.get_query_names()
                    seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                    output_seq[pcol.pos+1]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}
                    output_seq[pcol.pos+1]['ref']=base_to_num_map[r]
    
    return output_seq

def get_snp_training_pileups(dct):
    
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    
    include_intervals = None
    
    tbx = pysam.TabixFile(dct['include_bed'])
    include_intervals=IntervalTree(Interval(int(row[1]), int(row[2]), "%s" % (row[1])) for row in tbx.fetch(chrom, parser=pysam.asBed()))

    
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    vcf_path=dct['vcf_path']

    bcf_in = VariantFile(vcf_path)
    
    gt_map={(0,0):0, (1,1):0, (2,2):0, (1,2):1, (2,1):1, (0,1):1, (1,0):1, (0,2):1,(2,0):1}
    tr_pos={}
    for rec in bcf_in.fetch(chrom,start,end+1):
        gt=rec.samples.items()[0][1].get('GT')
        if gt in gt_map:
            if base_to_num_map[rec.alleles[gt[0]]]<4 and base_to_num_map[rec.alleles[gt[1]]]<4:
                tr_pos[rec.pos]=(gt_map[gt],base_to_num_map[rec.alleles[gt[0]]],base_to_num_map[rec.alleles[gt[1]]])
        
        
    nbr_size=20
    
    cnd_seq_freq=get_nbr(dct, nbr_type='freq')
    cnd_pos_freq=np.array(list(cnd_seq_freq.keys()))
    
    cnd_seq_gtruth=get_nbr(dct, nbr_type='gtruth')
    cnd_pos_gtruth=np.array(list(cnd_seq_gtruth.keys()))
    
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    
    pileup_dict={}
    
    output={'pos':[],0:[],5:[],10:[],15:[],20:[],25:[]}
    
    flag=0x4|0x100|0x200|0x400|0x800
        
    for pcol in samfile.pileup(chrom,max(0,start-1),end,min_base_quality=0,\
                                           flag_filter=flag,truncate=True):
            
            r=ref_dict[pcol.pos+1]
            if in_bed(include_intervals, pcol.pos+1) and r in 'AGTC':
                n=pcol.get_num_aligned()
                
                if n<=dct['maxcov'] and n>=dct['mincov'] and pcol.pos+1>=start and pcol.pos+1<=end:
                    
                    
                    if pcol.pos+1 in tr_pos:
                        seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                        name=pcol.get_query_names()
                        pileup_dict[pcol.pos+1]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}
                        output['pos'].append(pcol.pos+1)
                    
                    else:
                        seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                        alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=r and x[0] in 'AGTC')]+[0])/n
                        
                        if alt_freq>=0.10:
                            name=pcol.get_query_names()
                            pileup_dict[pcol.pos+1]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}

                            if 0.10<=alt_freq<0.15:
                                output[10].append(pcol.pos+1)
                            elif 0.15<=alt_freq<0.20:
                                output[15].append(pcol.pos+1)
                            elif 0.20<=alt_freq<0.25:
                                output[20].append(pcol.pos+1)
                            elif 0.25<=alt_freq:
                                output[25].append(pcol.pos+1)
                    
                        elif np.random.randint(2):
                            seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                            name=pcol.get_query_names()

                            pileup_dict[pcol.pos+1]={n:base_to_num_map[s] for (n,s) in zip(name,seq)}

                            if alt_freq<0.05:
                                output[0].append(pcol.pos+1)
                            
                            elif 0.05<=alt_freq<0.10:
                                output[5].append(pcol.pos+1)
                            
    
    pileup_list={'pos':[],'neg':[]}

    if output['pos']:
        tr_len=len(output['pos'])    
    else:
        tr_len=1e16
        
    sizes={0:tr_len, 5:tr_len//3,10:tr_len//3,15:tr_len//3, 20:tr_len, 25:tr_len}

    
    for instance in ['pos',0,5,10,15,20,25]:
        pos_list=output[instance]
        
        if pos_list:
            if instance!='pos':
                if sizes[instance]<len(output[instance]):
                    perm=np.random.permutation(sizes[instance])
                    pos_list=np.take(pos_list,perm,axis=0)
        
            for v_pos in pos_list:
                for nbr_types in [(cnd_pos_freq, cnd_seq_freq),(cnd_pos_gtruth, cnd_seq_gtruth)]:
                    cnd_pos, cnd_seq = nbr_types
                    ls_total_1, ls_total_2 = get_cnd_pos(v_pos, cnd_pos, dct['seq'])

                    nbr_dict={}

                    rlist=[(v_pos,base_to_num_map[ref_dict[v_pos]]) ]

                    sample=pileup_dict[v_pos].keys()

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
                                rlist.append((nb_pos, cnd_seq[nb_pos]['ref']))

                    rlist=sorted(rlist, key=lambda x:x[0])
                    total_rlist=np.array([x[1] for x in rlist])

                    mat=np.dstack([np.sum(np.eye(5)[tmp_mat[tmp_mat[:,len(ls_total_1)]==i]],axis=0) for i in range(4)]).transpose(2,0,1)[:,:,:4]

                    total_ref=np.eye(5)[total_rlist.astype(int)]
                    total_ref[:,4]=0
                    total_ref=total_ref[np.newaxis,:]
                    mat=np.dstack([mat,np.zeros([4,mat.shape[1]])+np.eye(4)[base_to_num_map[ref_dict[v_pos]]][:,np.newaxis]])
                    data=np.vstack([total_ref,np.multiply(mat,1-2*total_ref)])
                    data=np.hstack([np.zeros([5,nbr_size-len(ls_total_1),5]),data,np.zeros([5, nbr_size-len(ls_total_2), 5])]).astype(np.int8)

                    if instance=='pos':
                        pileup_list['pos'].append((v_pos,tr_pos[v_pos][0],tr_pos[v_pos][1],tr_pos[v_pos][2], base_to_num_map[ref_dict[v_pos]],data))

                    else:
                        pileup_list['neg'].append((v_pos,0,base_to_num_map[ref_dict[v_pos]], base_to_num_map[ref_dict[v_pos]], base_to_num_map[ref_dict[v_pos]],data))
               

    return pileup_list, dct['type']

def generate(params,mode='training'):
    cores=params['cpu']

    chrom_list=params['chrom_list']
    chrom_lengths=params['chrom_lengths']
    
    threshold=params['threshold']
    
    
    pool = mp.Pool(processes=args.cpu)
    
    in_dict_list=[]
    
    r_count=0
    for chrom in chrom_list:
        start, end=1, chrom_lengths[chrom]
        for mbase in range(start, end,1000000):
            r_count+=1
            d = copy.deepcopy(params)
            d['chrom']=chrom
            d['start']=mbase
            d['end']=min(end,mbase+1000000)
            d['type']='test' if r_count%40==0 else 'train'
            in_dict_list.append(d)
        
    results_dict=pool.imap_unordered(get_snp_training_pileups, in_dict_list)
    
    count=0
    file_count=0
    
    pileup_file_name=os.path.join(params['out_path'],'train.pileups.%d' %file_count) 
    pileup_file=open(pileup_file_name , "w")
    
    test_pileup_file_name=os.path.join(params['out_path'],'test.pileups') 
    test_pileup_file=open(test_pileup_file_name , "w")
    
    total=len(in_dict_list)
        
    for result_pair in results_dict:
        result, data_type=result_pair
        
        count+=1
        if count%100==0:
            pileup_file.close()
            file_count+=1

            pileup_file_name=os.path.join(params['out_path'],'train.pileups.%d' %file_count)
            pileup_file=open(pileup_file_name , "w")
            
        for i in ['pos','neg']:
            pileups=result[i]
            if pileups:
                if data_type=='train':
                    for data in pileups:
                            pos,gt,allele1,allele2,ref,mat=data
                            mat=mat.reshape(-1)
                            s='%s%d%d%d%d%d%s' %((11-len(str(pos)))*'0',pos,gt,allele1,allele2,ref,''.join([(6-len(x))*' '+x for x in mat.astype('<U6')]))
                            pileup_file.write(s)
                else:
                    for data in pileups:
                            pos,gt,allele1,allele2,ref,mat=data
                            mat=mat.reshape(-1)
                            s='%s%d%d%d%d%d%s' %((11-len(str(pos)))*'0',pos,gt,allele1,allele2,ref,''.join([(6-len(x))*' '+x for x in mat.astype('<U6')]))
                            test_pileup_file.write(s)
        print('%s: %d/%d regions done\n'%(str(datetime.datetime.now()),count,total),flush=True)

    results_dict=None
    elapsed=time.time()-t
    
    pileup_file.close()
    test_pileup_file.close()
    
if __name__ == '__main__':
    
    t=time.time()

    print('%s: Starting pileup generation.' %str(datetime.datetime.now()))
    parser = argparse.ArgumentParser()

    #-r chromosome region   -m mode   -bam bam file   -ref reference file   -vcf ground truth variants   -o output path
    
    parser.add_argument("-chrom", "--chrom", help="Chromosome region")
    parser.add_argument("-seq",  "--sequencing",  help="Sequencing type, options are 'ont' and 'pacbio'", type=str, default='ont')
    
    parser.add_argument("-bam", "--bam", help="Bam file")
    parser.add_argument("-ref", "--ref", help="Size")
    parser.add_argument("-vcf", "--vcf", help="Ground truth variants")
    
    parser.add_argument("-o", "--output", help="Output path")
    
    parser.add_argument("-cpu", "--cpu", help="CPUs",type=int)
    parser.add_argument("-include_bed", "--include_bed", help="Include BED file")
    parser.add_argument("-exclude_bed", "--exclude_bed", help="Exclude BED file")

    parser.add_argument("-start", "--start", help="start",type=int)
    parser.add_argument("-end", "--end", help="end",type=int)
    
    parser.add_argument("-nbr_t",  "--neighbor_threshold",  help="SNP neighboring site thresholds with lower and upper bounds seperated by comma, for Nanopore reads '0.4,0.6' is recommended, for PacBio CCS anc CLR reads '0.3,0.7' and '0.3,0.6' are recommended respectively", type=str, default='0.4,0.6')
    
    parser.add_argument("-mincov",  "--mincov",  help="min coverage", type=int, default=8)
    parser.add_argument("-maxcov",  "--maxcov",  help="max coverage", type=int, default=160)
    
    parser.add_argument('-wgs_contigs_type','--wgs_contigs_type', \
                        help="""Options are "with_chr", "without_chr" and "all",\
                        "with_chr" option will assume \
                        human genome and run NanoCaller on chr1-22, "without_chr" will \
                        run on chromosomes 1-22 if the BAM and reference genome files \
                        use chromosome names without "chr". "all" option will run \
                        NanoCaller on each contig present in reference genome FASTA file.""", \
                        type=str, default='with_chr')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    if args.chrom:
        chrom_list= args.chrom.split()
    
    else:
        if args.wgs_contigs_type=='with_chr':
            chrom_list=['chr%d' %d for d in range(1,23)]

        elif args.wgs_contigs_type == 'without_chr':
            chrom_list=['%d' %d for d in range(1,23)]

        elif args.wgs_contigs_type == 'all':
            chrom_list=[]

            try:
                with open(args.ref+'.fai','r') as file:
                    for line in file:
                        chrom_list.append(line.split('\t')[0])

            except FileNotFoundError:
                print('%s: index file .fai required for reference genome file.\n' %str(datetime.datetime.now()), flush=True)
                sys.exit(2)
    
    chrom_lengths={}
    with open(args.ref+'.fai','r') as file:
        for line in file:
            chrom_lengths[line.split('\t')[0]]=int(line.split('\t')[1])
            
            
    threshold=[float(args.neighbor_threshold.split(',')[0]), float(args.neighbor_threshold.split(',')[1])]
    
    
    in_dict={'chrom_list':chrom_list,'chrom_lengths':chrom_lengths,\
         'sam_path':args.bam, 'fasta_path':args.ref, 'vcf_path':args.vcf,\
             'out_path':args.output,'seq':args.sequencing, \
             'threshold':threshold, 'cpu':args.cpu, 'include_bed':args.include_bed,'exclude_bed':args.exclude_bed,\
            'mincov':args.mincov,'maxcov':args.maxcov}    
    
    generate(in_dict)
    elapsed=time.time()-t
    print ('%s: Total Time Elapsed: %.2f seconds' %(str(datetime.datetime.now()), elapsed))