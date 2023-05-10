import sys,pysam, time,os,copy,argparse,subprocess, random, re, datetime
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from subprocess import Popen, PIPE, STDOUT
from Bio import pairwise2
from intervaltree import Interval, IntervalTree

gt_map={(0,0):0, (1,1):1, (0,1):2, (1,0):2,(1,2):3, (2,1):3}

mapping={'A':0,'G':1,'T':2,'C':3,'-':4}
rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

allele_map={('I','I'):0, ('D','D'):1, ('N','I'):2, ('I','N'):3, ('N','D'):4, ('D','N'):5, \
           ('I','D'):6, ('D','I'):7, ('N','N'):8}

def pairwise(x,y):
    alignments = pairwise2.align.globalms(x, y, 2, -1.0, -0.9, -0.1)

    return alignments

def msa(seq_list, ref, v_pos, mincov, maxcov):
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
        return (0,0,None,None,None)
        
    try:
        ref_real_0_mat=np.eye(5)[[mapping[x] for x in ref_real_0[:128]]].transpose()
    except UnboundLocalError:
        return (0,0,None)
    
    mat=np.array([[mapping[c] for c in x] for x in zz_0])
    h0_mat=np.sum(np.eye(5)[mat],axis=0).transpose()

    return (1,1,np.dstack([h0_mat, ref_real_0_mat]))

def get_indel_training_pileups(dct):
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
    mincov,maxcov=dct['mincov'],dct['maxcov']
    
    include_intervals, exclude_intervals=None, None
    
    
    if dct['include_bed']:
        tbx = pysam.TabixFile(dct['include_bed'])
        include_intervals=IntervalTree(Interval(int(row[1]), int(row[2]), "%s" % (row[1])) for row in tbx.fetch(chrom, parser=pysam.asBed()))
        
        def in_bed(tree,pos):            
            return tree.overlaps(pos)
        
        include_intervals=IntervalTree(include_intervals.overlap(start,end))
        
        if not include_intervals:
            return None
    
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
    
    bcf_in = VariantFile(dct['vcf_path'])
    tr_var={}
    for rec in bcf_in.fetch(chrom,start,end+1):
            gt=rec.samples.items()[0][1].get('GT')
            tr_var[rec.pos-10]=gt_map[gt]
            tr_var[rec.pos]=gt_map[gt]
            
    tr_pos=set(tr_var.keys())        

    output={'pos':[],'high':[],'low':[]}

    for pcol in samfile.pileup(chrom,max(0,start-100),end,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
        v_pos=pcol.pos+1
        if v_pos in tr_pos:
            read_names=pcol.get_query_names()
            read_names_0=set(read_names) & hap_reads_0
            read_names_1=set(read_names) & hap_reads_1
            len_seq_0=len(read_names_0)
            len_seq_1=len(read_names_1)
            
            if len(read_names_0)>=dct['mincov'] and len(read_names_1)>=dct['mincov']:
                output['pos'].append(pcol.pos+1)
                        
        if in_bed(include_intervals, v_pos) and not ex_bed(exclude_intervals, v_pos):                
            read_names=pcol.get_query_names()
            read_names_0=set(read_names) & hap_reads_0
            read_names_1=set(read_names) & hap_reads_1
            len_seq_0=len(read_names_0)
            len_seq_1=len(read_names_1)
            
            if len_seq_0>=dct['mincov']  and len_seq_1>=dct['mincov']:
                seq=[x for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)]

                tmp_seq_0=''.join([s for n,s in zip(read_names,seq) if n in read_names_0])
                tmp_seq_1=''.join([s for n,s in zip(read_names,seq) if n in read_names_1])

                del_freq_0=(tmp_seq_0.count('-'))/len_seq_0 if len_seq_0>0 else 0
                ins_freq_0=tmp_seq_0.count('+')/len_seq_0 if len_seq_0>0 else 0

                del_freq_1=(tmp_seq_1.count('-'))/len_seq_1 if len_seq_1>0 else 0
                ins_freq_1=tmp_seq_1.count('+')/len_seq_1 if len_seq_1>0 else 0

                if 0.3<=del_freq_0 or 0.3<=del_freq_1 or 0.3<=ins_freq_0 or 0.3<=ins_freq_1:
                    output['high'].append(pcol.pos+1)

                elif del_freq_0<=0.2 and del_freq_1<=0.2 and ins_freq_0<=0.1 and ins_freq_1<=0.1 and np.random.randint(100)==0:
                    output['low'].append(pcol.pos+1)
    
    
    if output['pos']:
        tr_len=len(output['pos'])    
    else:
        tr_len=20
        
    sizes={'high':tr_len, 'low':tr_len}

    output['pos']=set(output['pos'])
    
    
    for i in ['high','low']:
        if sizes[i]<len(output[i]):
            perm=np.random.permutation(sizes[i])
            var_list=np.take(output[i],perm,axis=0)
            output[i]=set(var_list)
            
            
    
    total_variants=set.union(output['pos'],output['high'],output['low'])
    
    data=[]

    for pcol in samfile.pileup(chrom,max(0,start-100),end,min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
        
        v_pos=pcol.pos+1
        
        if v_pos in total_variants:
            read_names=pcol.get_query_names()
            read_names_0=set(read_names) & hap_reads_0
            read_names_1=set(read_names) & hap_reads_1
           

            d={'hap0':{},'hap1':{}}
            d_tot={}
                        
            ref=''.join([ref_dict[p] for p in range(v_pos-window_before,v_pos+window_after+1)])
            
            if len(ref)<window_after:
                continue
                    
            for pread in pcol.pileups:
                dt=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]                    
                d_tot[pread.alignment.qname]=dt
                    
                if pread.alignment.qname in read_names_0:
                    d['hap0'][pread.alignment.qname]=dt

                elif pread.alignment.qname in read_names_1:
                    d['hap1'][pread.alignment.qname]=dt    
            
            seq_list=d['hap0']
            flag0,indel_flag0,data_0=msa(seq_list,ref,v_pos,2,dct['maxcov'])

            seq_list=d['hap1']
            flag1,indel_flag1,data_1=msa(seq_list,ref,v_pos,2,dct['maxcov'])
            
            seq_list = d_tot
            flag_total,indel_flag_total,data_total=msa(seq_list,ref,v_pos,dct['mincov'],dct['maxcov'])

            if flag0 and flag1 and flag_total:
                gt= tr_var[v_pos] if v_pos in tr_pos else 0
                data.append((v_pos,gt,(data_0,data_1,data_total)))
    
    return data

def generate(params,mode='training'):
    cores=params['cpu']

    chrom_list=params['chrom_list']
    chrom_lengths=params['chrom_lengths']    
    
    pool = mp.Pool(processes=args.cpu)
    
    in_dict_list=[]
    
    r_count=0

    for chrom in chrom_list:
        start, end=1, chrom_lengths[chrom]
        for mbase in range(start, end,100000):
            d = copy.deepcopy(params)
            d['chrom']=chrom
            d['start']=mbase
            d['end']=min(end,mbase+100000)
            in_dict_list.append(d)
            
    results_dict=pool.imap_unordered(get_indel_training_pileups, in_dict_list)
    
    count=0
    file_count=0
    
    pileup_file_name=os.path.join(params['out_path'],'train.pileups.%d' %file_count) 
    pileup_file=open(pileup_file_name , "w")
    
    test_pileup_file_name=os.path.join(params['out_path'],'test.pileups') 
    test_pileup_file=open(test_pileup_file_name , "w")
    
    total=len(in_dict_list)
    for result in results_dict:        
        count+=1
        if count%1000==0:
            pileup_file.close()
            file_count+=1

            pileup_file_name=os.path.join(params['out_path'],'train.pileups.%d' %file_count)
            pileup_file=open(pileup_file_name , "w")
            
        if result:
            if count%40!=0:
                for data in result:
                    pos,gt,phased_data=data
                    mat_0=phased_data[0].reshape(-1).astype(np.int16)
                    mat_1=phased_data[1].reshape(-1).astype(np.int16)
                    mat_2=phased_data[2].reshape(-1).astype(np.int16)
                    s='%s%d%d%s%s%s' %((11-len(str(pos)))*'0',pos,gt, ''.join([(4-len(x))*' '+x for x in mat_0.astype('<U4')]), ''.join([(4-len(x))*' '+x for x in mat_1.astype('<U4')]),''.join([(4-len(x))*' '+x for x in mat_2.astype('<U4')]))
                    pileup_file.write(s)
            else:
                for data in result:
                    pos,gt,phased_data=data
                    mat_0=phased_data[0].reshape(-1).astype(np.int16)
                    mat_1=phased_data[1].reshape(-1).astype(np.int16)
                    mat_2=phased_data[2].reshape(-1).astype(np.int16)
                    s='%s%d%d%s%s%s' %((11-len(str(pos)))*'0',pos,gt, ''.join([(4-len(x))*' '+x for x in mat_0.astype('<U4')]), ''.join([(4-len(x))*' '+x for x in mat_1.astype('<U4')]),''.join([(4-len(x))*' '+x for x in mat_2.astype('<U4')]))
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
