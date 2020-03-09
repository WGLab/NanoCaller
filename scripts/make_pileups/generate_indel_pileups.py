import sys,pysam, time,os,re,copy,argparse,gzip,itertools,subprocess,git
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree
from Bio import SeqIO
from Bio.Align.Applications import TCoffeeCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from tempfile import mkstemp
from subprocess import Popen, PIPE, STDOUT

os.environ["PATH"] += os.pathsep + '/home/ahsanm/lib/tcoffee/t_coffee/src'
window_before,window_after=40,100
gt_map={(0,0):0, (1,1):0, (1,2):2, (2,1):2, (0,1):1, (1,0):1}

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

    '''filter_process = Popen(['t_coffee','-other_pg','seq_reformat', '-in','stdin','-action','+trim','_seq_O85','+upper','-output' ,'fasta'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    fltr_seqs=filter_process.communicate(input=fa_tmp_file.encode('utf-8'))[0]

    fltr_seqs+=b'>ref_SEQ\n'
    fltr_seqs+= b'%b\n' %bytes(ref.encode('utf-8'))'''

    fa_tmp_file+='>ref_SEQ\n'
    fa_tmp_file+= '%s' %ref
    
    msa_process = Popen(['muscle', '-quiet'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    hap_file=msa_process.communicate(input=fa_tmp_file.encode('utf-8'))

    if len(hap_file)==0:
        print('hapfile length 0')


    tmp=hap_file[0].decode('utf-8')[1:].replace('\n','').split('>')

    zz_0=[]
    for seq in tmp:
        p1,p2=seq.split('_SEQ')
        if p1!='ref':
            zz_0.append(p2)
        else:
            ref_real_0=p2

    if len(zz_0)<mincov:
        return (0,0,None)
    cnt=0
    cnt_del=0
    seq_ts=ref_real_0

    while cnt<window_before or seq_ts[cnt_del]=='-':
        if seq_ts[cnt_del]!='-':
            cnt+=1
        cnt_del+=1

    ref_real_0_mat=np.eye(5)[[mapping[x] for x in ref_real_0[cnt_del-1:cnt_del+63]]].transpose()
    
    mat=np.array([[mapping[c] for c in x] for x in zz_0])
    try:
        h0_mat=np.sum(np.eye(5)[mat],axis=0).transpose()
    except IndexError:
        print(v_pos)
        print('ERROR FOUND @@@@@@@@@@@@@@@@@@@@@@@@@@@@|||||||||@@@@@@@@@@@@@@@@@@@@')
        print(fa_tmp_file)
        print('-------------------SEQ LISTTTTT------------')
        print(seq_list)
        print(mat)
        return
    
    h0_mat=h0_mat[:,cnt_del-1:cnt_del+63]
    h0_mat=h0_mat
    
    h0_mat_tmp=h0_mat.astype(np.float32)
    h0_mat_tmp=h0_mat_tmp/(np.sum(h0_mat_tmp,axis=0))-ref_real_0_mat
    
    indel_flag=np.sum(np.abs(h0_mat_tmp[4,:])>=0.5)>0
    
    #final_mat_0=ref_real_0_mat-h0_mat
    
    return (1,indel_flag,np.dstack([h0_mat, ref_real_0_mat])) 
    

def get_training_candidates(dct):
    chrom=dct['chrom']
    start=dct['start']
    end=dct['end']
    sam_path=dct['sam_path']
    fasta_path=dct['fasta_path']
    vcf_path=dct['vcf_path']
    bed_path=dct['bed']
        
    if bed_path:
        bed_file=os.path.join(bed_path,'%s.bed' %chrom)
        with open(bed_file) as file:
            content=[x.rstrip('\n') for x in file]

        content=[x.split('\t')[1:] for x in content]
        content=[(int(x[0]),int(x[1])) for x in content]
        t=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)


    samfile = pysam.Samfile(sam_path, "rb")
    hap_dict={1:[],2:[]}
    for pread in samfile.fetch(chrom, dct['start']-1, dct['end']+1):
        if pread.has_tag('HP'):
            hap_dict[pread.get_tag('HP')].append(pread.qname)

    hap_reads_0=set(hap_dict[1])
    hap_reads_1=set(hap_dict[2])
    

    bcf_in = VariantFile(dct['vcf_path'])
    tr_pos={}
    for rec in bcf_in.fetch(chrom,start,end+1):
        if rec.samples.items()[0][1].get('PS') or rec.samples.items()[0][1].get('GT') in [(1,1),(1,2),(2,1)]:
            gt=rec.samples.items()[0][1].get('GT')
            vt=v_type(rec.ref,rec.alleles[gt[0]],rec.alleles[gt[1]])
            if vt!=8:
                tr_pos[rec.pos]=vt

        
    fastafile=pysam.FastaFile(fasta_path)
    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1),\
                                                                 fastafile.fetch(chrom,max(1,start-40)-1,end+40)) }
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
    
    pileup_dict={}
    output={'pos':[],'high':[],'low':[]}

    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            
            make=False
            if ref_df.loc[pcol.pos+1].ref !='*':
                                    
                if bed_path:
                    if not t[pcol.pos+1]:
                        continue
                        
                if pcol.pos+1 in tr_pos.keys():
                    read_names=pcol.get_query_names()
                    read_names_0=[x for x in read_names if x in hap_reads_0]
                    read_names_1=[x for x in read_names if x in hap_reads_1]
                    
                    if len(read_names_0)>=dct['mincov'] and len(read_names_1)>=dct['mincov']:
                        make=True
                        output['pos'].append(pcol.pos+1,)
                    
                elif np.random.randint(20)==0:
                    read_names=pcol.get_query_names()
                    read_names_0=[x for x in read_names if x in hap_reads_0]
                    read_names_1=[x for x in read_names if x in hap_reads_1]
                    
                    if len(read_names_0)>=dct['mincov'] and len(read_names_1)>=dct['mincov']:
                        seq=[x[:2].upper() for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)]
                    
                        tmp_0=[s for n,s in zip(read_names,seq) if n in read_names_0]
                        len_seq_0=len(tmp_0)
                        tmp_seq_0=''.join(tmp_0)

                        tmp_1=[s for n,s in zip(read_names,seq) if n in read_names_1]
                        len_seq_1=len(tmp_1)
                        tmp_seq_1=''.join(tmp_1)

                        del_freq_0=tmp_seq_0.count('-')/len_seq_0 if len_seq_0>0 else 0
                        ins_freq_0=tmp_seq_0.count('+')/len_seq_0 if len_seq_0>0 else 0

                        del_freq_1=tmp_seq_1.count('-')/len_seq_1 if len_seq_1>0 else 0
                        ins_freq_1=tmp_seq_1.count('+')/len_seq_1 if len_seq_1>0 else 0


                        if (0.6<=del_freq_0 or 0.6<=del_freq_1 or 0.3<=ins_freq_0 or 0.3<=ins_freq_1):
                            make=True
                            
                            output['high'].append(pcol.pos+1)

                        elif (del_freq_0<=0.2 and del_freq_1<=0.2 and ins_freq_0<=0.1 and ins_freq_1<=0.1) and np.random.randint(100)==0:
                            make=True
                            output['low'].append(pcol.pos+1)
                if make:
                    d={'hap0':{},'hap1':{}}
                    for pread in pcol.pileups:
                        if pread.alignment.qname in read_names_0:
                            d['hap0'][pread.alignment.qname]=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]

                        elif pread.alignment.qname in read_names_1:
                            d['hap1'][pread.alignment.qname]=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]
                        
                    pileup_dict[pcol.pos+1]=d
    pileup_list={'pos':[],'high':[],'low':[]}
    
    np.random.seed(76)
    
    if output['pos']:
        tr_len=len(output['pos'])    
    else:
        tr_len=20
        
    sizes={'high':tr_len, 'low':tr_len}

    for i in ['pos','high','low']:
        pos_list=output[i]
        if pos_list:
            if i!='pos':
                if sizes[i]<len(output[i]):
                    perm=np.random.permutation(sizes[i])
                    pos_list=np.take(pos_list,perm,axis=0)
            
            for v_pos in pos_list:
                d=pileup_dict[v_pos]
                ref=''.join(ref_df.reindex(range(v_pos-window_before,v_pos+window_after+1)).fillna('').ref.to_list())
                seq_list=d['hap0']
                flag0,indel_flag0,data_0=msa(seq_list,ref,v_pos,dct['mincov'])

                seq_list=d['hap1']
                flag1,indel_flag1,data_1=msa(seq_list,ref,v_pos,dct['mincov'])
                
                if flag0 and flag1:
                    data=(data_0,data_1)
                    if i=='pos':
                        pileup_list[i].append((v_pos,tr_pos[v_pos],data))
                    else:
                        pileup_list[i].append((v_pos,8,data))
                
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
    
    hap_dict={1:[],2:[]}
    
    for pread in samfile.fetch(chrom, max(0,dct['start']-100), dct['end']+100):
        if pread.has_tag('HP'):
            hap_dict[pread.get_tag('HP')].append(pread.qname)

    hap_reads_0=set(hap_dict[1])
    hap_reads_1=set(hap_dict[2])
    
    
    pileup_list=[]
    
    prev=0
    for pcol in samfile.pileup(chrom,max(0,start-1),end,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            
            
            
            if pcol.pos+1>prev+10:                
                               

                if bed_path:
                    if not t[pcol.pos+1]:
                        continue
                        
                        
                read_names=pcol.get_query_names()
                read_names_0=[x for x in read_names if x in hap_reads_0]
                read_names_1=[x for x in read_names if x in hap_reads_1]

                if len(read_names_0)>=dct['mincov'] and len(read_names_1)>=dct['mincov']:
                    seq=[x[:2].upper() for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False, add_indels=True)]
                    
                    tmp_0=[s for n,s in zip(read_names,seq) if n in read_names_0]
                    len_seq_0=len(tmp_0)
                    tmp_seq_0=''.join(tmp_0)
                    
                    tmp_1=[s for n,s in zip(read_names,seq) if n in read_names_1]
                    len_seq_1=len(tmp_1)
                    tmp_seq_1=''.join(tmp_1)
                                        
                    del_freq_0=tmp_seq_0.count('-')/len_seq_0 if len_seq_0>0 else 0
                    ins_freq_0=tmp_seq_0.count('+')/len_seq_0 if len_seq_0>0 else 0
                    
                    del_freq_1=tmp_seq_1.count('-')/len_seq_1 if len_seq_1>0 else 0
                    ins_freq_1=tmp_seq_1.count('+')/len_seq_1 if len_seq_1>0 else 0
                    
                    
                    if len_seq_0>=dct['mincov']  and len_seq_1>=dct['mincov']  and (0.4<=del_freq_0 or 0.4<=del_freq_1 or 0.25<=ins_freq_0 or 0.25<=ins_freq_1):
                        v_pos=pcol.pos+1
                        d={'hap0':{},'hap1':{}}
                        for pread in pcol.pileups:
                            if pread.alignment.qname in read_names_0:
                                d['hap0'][pread.alignment.qname]=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]

                            elif pread.alignment.qname in read_names_1:
                                d['hap1'][pread.alignment.qname]=pread.alignment.query_sequence[max(0,pread.query_position_or_next-window_before):pread.query_position_or_next+window_after]
                                
                        ref=''.join(ref_df.reindex(range(v_pos-window_before,v_pos+window_after+1)).fillna('').ref.to_list())
                        seq_list=d['hap0']
                        flag0, indel_flag0, data_0=msa(seq_list,ref,v_pos,dct['mincov'])

                        
                        seq_list=d['hap1']
                        flag1,indel_flag1,data_1=msa(seq_list,ref,v_pos,dct['mincov'])

                        if flag0 and flag1  and (indel_flag0 or indel_flag1):
                            data=(data_0,data_1)
                            prev=pcol.pos+1
                            
                            pileup_list.append((v_pos,data))
                                
    return pileup_list
    
    

def generate(params,mode='training'):
    cores=params['cpu']
    mode=params['mode']
    chrom=params['chrom']
    threshold=params['threshold']
    start,end=params['start'],params['end']
    
    if mode in ['training','train']: #'pos,ref,seq,names'
        print('starting training pileups',flush=True)

        pool = mp.Pool(processes=cores)
        fname='%s.pileups.' %params['chrom']
        file_list={}
        for i in ['pos','high','low']:
            suffix='pos' if i == 'pos' else 'neg.%s' %i
            tmp_name=os.path.join(params['out_path'],'%s%s' %(fname,suffix)) 
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
                #(v_pos,ref_df.loc[v_pos].ref,data)
                #(v_pos,tr_pos[v_pos],data)    output=(gt,ref[0],allele_1,allele_2)
            for result in results_dict:
                for i in ['pos','high','low']:
                    pileups=result[i]
                    if pileups:
                        
                        for data in pileups:
                                pos,gt,phased_data=data
                                mat_0=phased_data[0].reshape(-1).astype(np.int16)
                                mat_1=phased_data[1].reshape(-1).astype(np.int16)

                                s='%s%d%d%s%s' %((11-len(str(pos)))*'0',pos,gt, ''.join([(4-len(x))*' '+x for x in mat_0.astype('<U4')]), ''.join([(4-len(x))*' '+x for x in mat_1.astype('<U4')]))
                                file_list[i].write(s)
            
            results_dict=None
            elapsed=time.time()-t
            
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)                    

                
    elif mode in ['testing','test']:#'pos,depth,freq,ref,seq,names'
        
        print('starting pileups',flush=True)
        pool = mp.Pool(processes=cores)
        
        fname='%s.pileups.test' %params['chrom']
        fname=os.path.join(params['out_path'],fname)
        file=open(fname , "w")
        
        start,end=params['start'],params['end']
        
        for mbase in range(start,end,int(1e7)):
            print('starting pool:'+str(mbase),flush=True)
            t=time.time()                

            in_dict_list=[]
            for k in range(mbase,min(end,mbase+int(1e7)),10000):
                d = copy.deepcopy(params)
                d['start']=k
                d['end']=min(end,k+10000)
                in_dict_list.append(d)
            results_dict = pool.map(get_testing_candidates, in_dict_list)
            
            for result in results_dict:
                    if result:
                        
                        for data in result:
                            pos,phased_data=data
                            mat_0=phased_data[0].reshape(-1).astype(np.int16)
                            mat_1=phased_data[1].reshape(-1).astype(np.int16)

                            s='%s%d%s%s' %((11-len(str(pos)))*'0',pos, ''.join([(4-len(x))*' '+x for x in mat_0.astype('<U4')]), ''.join([(4-len(x))*' '+x for x in mat_1.astype('<U4')]))
                            
                            file.write(s)

            results_dict=None
            elapsed=time.time()-t
            print ('Elapsed: %.2f seconds' %elapsed,flush=True)
            print('finishing pool:'+str(mbase),flush=True)
            
    
if __name__ == '__main__':
    print('git commit hash: %s' %str(git.Repo("/home/ahsanm/repos/NanoVar").heads[0].commit))
    chrom_length={'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555, 'chr5':181538259, 'chr6':170805979, \
             'chr7':159345973, 'chr8':145138636, 'chr9':138394717, 'chr10':133797422, 'chr11':135086622, 'chr12':133275309,\
             'chr13':114364328, 'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,\
             'chr19':58617616, 'chr20':64444167, 'chr21':46709983, 'chr22':50818468, 'chrX':156040895, 'chrY':57227415,'1':248956422, '2':242193529, '3':198295559, '4':190214555, '5':181538259, '6':170805979,'7':159345973, '8':145138636, '9':138394717, '10':133797422, '11':135086622, '12':133275309,'13':114364328, '14':107043718, '15':101991189, '16':90338345, '17':83257441, '18':80373285,'19':58617616, '20':64444167, '21':46709983, '22':50818468, 'X':156040895, 'Y':57227415}
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
    parser.add_argument("-test_vcf", "--test_vcf", help="Test VCF")
    parser.add_argument("-nbr_vcf", "--nbr_vcf", help="VCF for generating neighbors")

    
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
            'mincov':args.mincov, 'nbr_size':args.nbr_size, 'test_vcf':args.test_vcf, 'nbr_vcf': args.nbr_vcf}    
    
    t=time.time()
    generate(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)