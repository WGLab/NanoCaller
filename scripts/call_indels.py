import sys,pysam, time,os,re,copy,argparse,gzip,itertools,subprocess
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
chrom_length={'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555, 'chr5':181538259, 'chr6':170805979, \
             'chr7':159345973, 'chr8':145138636, 'chr9':138394717, 'chr10':133797422, 'chr11':135086622, 'chr12':133275309,\
             'chr13':114364328, 'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,\
             'chr19':58617616, 'chr20':64444167, 'chr21':46709983, 'chr22':50818468, 'chrX':156040895, 'chrY':57227415}

mapping={'A':0,'G':1,'T':2,'C':3,'*':4,'N':4}

def pairwise(x,y):
    alignments = pairwise2.align.localms(x, y, 2, -1.0, -0.5, -0.1)

    # Use format_alignment method to format the alignments in the list
    #center=''
    
    #print(format_alignment(*alignments[0]))

    return alignments


def repeat_check(x):
    if len(x)==0:
        return False
    count=1
    letter=x[0]
    i=1
    while count<4 and i<len(x):
        
        if x[i]==letter:
            count+=1
        else:
            letter=x[i]
            count=1
        i+=1
    return count>=4


def indel_calls(params):
    window_before,window_after=30,35
    fasta_path=params['fasta_path']
    fastafile=pysam.FastaFile(fasta_path)
    sam_path=params['sam_path']
    samfile = pysam.Samfile(sam_path, "rb")
    chrom=params['chrom']
    bed_path=params['bed']

    file=open(os.path.join(params['tree_path'],'%s.tree' %chrom),'r')
    content=file.readlines()
    file.close()
    
    hap_tree=IntervalTree()

    for line in content:
        line=line.rstrip('\n')
        line=line.split('>')
        begin,end=int(line[0].split('-')[0]),int(line[0].split('-')[1])
        hap_tree[begin:end+1]=set(line[1].split(':'))

    if bed_path:
        bed_file=os.path.join(bed_path,'%s.bed' %params['chrom'])
        
        file=open(bed_file)
        content=[x.rstrip('\n') for x in file]
        file.close()
        
        content=[x.split('\t')[1:] for x in content]
        content=[(int(x[0]),int(x[1])) for x in content]
        bed_tree=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)

    reads_dict={}
    for pread in samfile.fetch(chrom,params['start']-1,params['end']):
            reads_dict[pread.query_name]=pread.query_sequence

    output_list=[]

    prev_cand_reads={}        

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,params['start']-300),params['end']+300+1), fastafile.fetch(params['chrom'], max(1,params['start']-300)-1,params['end']+300)) }

    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)


    for pcol in samfile.pileup(chrom,params['start']-1,params['end'],min_base_quality=0,flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
        v_pos=pcol.pos+1
        if v_pos>=params['end']:
            pass#return output_list
        if bed_path:
            if not bed_tree[pcol.pos+1]:
                continue
        try:
            n=pcol.get_num_aligned()
            if n>params['mincov']:
                nbr=''.join(ref_df.reindex(range(pcol.pos-3,pcol.pos+6)).ref.to_list())
                check=repeat_check(nbr)

                if check==False:                    
                    seq=''.join(pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True))
                    p_count=seq.count('+')
                    n_count=seq.count('-')+seq.count('*')

                    if p_count/n>=0.2 or n_count/n>=0.6:
                        cand_reads={pread.alignment.query_name:pread.query_position for pread in pcol.pileups}
                        read_names=set(cand_reads.keys())

                        common_reads=[iv.data&read_names for iv in hap_tree[v_pos]]
                        common_reads=sorted(common_reads, key=lambda x: len(x),reverse=True)
                        if len(common_reads[:2])<2:
                            continue
                        best_1,best_2=common_reads[:2]

                        hap1_reads=[]
                        for tmp_read in best_1:
                            if reads_dict[tmp_read]:
                                if cand_reads[tmp_read]:
                                    rd=reads_dict[tmp_read][max(0,cand_reads[tmp_read]-window_before):min(cand_reads[tmp_read]+window_after+1,len(reads_dict[tmp_read])-1)]
                                    if rd:
                                        if len(rd)>4:
                                            hap1_reads.append(rd)
                                elif tmp_read in prev_cand_reads.keys():
                                    rd=reads_dict[tmp_read][max(0,prev_cand_reads[tmp_read]-window_before):min(prev_cand_reads[tmp_read]+window_after+1,len(reads_dict[tmp_read])-1)]
                                    if rd:
                                        if len(rd)>4:
                                            hap1_reads.append(rd)
                        hap2_reads=[]
                        for tmp_read in best_2:
                            if reads_dict[tmp_read]:
                                if cand_reads[tmp_read]:
                                    rd=reads_dict[tmp_read][max(0,cand_reads[tmp_read]-window_before):min(cand_reads[tmp_read]+window_after+1,len(reads_dict[tmp_read])-1)]
                                    if rd:
                                        if len(rd)>4:
                                            hap2_reads.append(rd)
                                elif tmp_read in prev_cand_reads.keys():
                                    rd=reads_dict[tmp_read][max(0,prev_cand_reads[tmp_read]-window_before):min(prev_cand_reads[tmp_read]+window_after+1,len(reads_dict[tmp_read])-1)]
                                    if rd:
                                        if len(rd)>4:
                                            hap2_reads.append(rd)


                        if len(hap1_reads)==0 or len(hap2_reads)==0:
                            continue
                        
                        hap1_dnd_tmp_file,hap1_dnd_tmp= mkstemp()
                        hap1_fa_tmp_file=''
                        seq_count=0
                        for b in hap1_reads:
                            hap1_fa_tmp_file+='>hap1seq%d\n'%seq_count
                            hap1_fa_tmp_file+= '%s\n' %b
                            seq_count+=1
                            
                        
                        hap2_dnd_tmp_file,hap2_dnd_tmp= mkstemp()
                        hap2_fa_tmp_file=''

                        seq_count=0
                        for b in hap2_reads:
                            hap2_fa_tmp_file+='>hap2seq%d\n'%seq_count
                            hap2_fa_tmp_file+='%s\n' %b
                            seq_count+=1

                        msa_process = Popen(['t_coffee', '-infile=stdin','-outfile','stdout','-out_lib=/dev/null', '-no_warning', '-quiet','-newtree=%s'%hap1_dnd_tmp,'-quicktree','-distance_matrix_mode=fast'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
                        hap1_file=msa_process.communicate(input=hap1_fa_tmp_file.encode('utf-8'))[0]
                        
                        try:
                            os.remove(hap1_dnd_tmp)
                        except FileNotFoundError:
                            pass
                        
                        hap1_re_al={}
                        hap1_file=hap1_file.decode('utf-8')
                        if 'hap' in hap1_file:
                            for line in hap1_file.split('\n'):
                                if 'hap' in line:
                                    if line.split()[0] not in hap1_re_al.keys():
                                        hap1_re_al[line.split()[0]]=''
                                    hap1_re_al[line.split()[0]]+=line.split()[1]

                        else:
                            continue
                            
                        seq_cons_nosub=''
                        seq_cons=''
                        seq_cons_nodel=''

                        seq_len=len(list(hap1_re_al.values())[0])

                        for i in range(seq_len):
                            counts={'A':0,'G':0,'T':0,'C':0,'-':0}
                            for read in hap1_re_al.values():
                                counts[read[i]]+=1
                            top_list=sorted(list(counts.items()), key=lambda x: x[1])
                            top=top_list[-1][0]
                            seq_cons_nosub+=top
                            if top=='-':
                                    if top_list[-2][1]/len(hap1_re_al.values())>=0.2:
                                        #print(top_list[-1][1],top_list[-2][1])
                                        top=top_list[-2][0]
                                        #print('subs %s for -' %top)
                                        #print(top_list)

                            seq_cons+=top
                            top ='' if top=='-' else top
                            seq_cons_nodel+=top

                        hap1_cons=seq_cons_nodel

                        msa_process = Popen(['t_coffee', '-infile=stdin','-outfile','stdout','-out_lib=/dev/null', '-no_warning', '-quiet','-newtree=%s'%hap1_dnd_tmp,'-quicktree','-distance_matrix_mode=fast'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
                        hap2_file=msa_process.communicate(input=hap2_fa_tmp_file.encode('utf-8'))[0]
                        
                        try:
                            os.remove(hap2_dnd_tmp)
                        except FileNotFoundError:
                            pass
                        
                        hap2_re_al={}
                        hap2_file=hap2_file.decode('utf-8')
                        if 'hap' in hap2_file:
                            for line in hap2_file.split('\n'):
                                if 'hap' in line or 'Ref' in line:
                                    if line.split()[0] not in hap2_re_al.keys():
                                        hap2_re_al[line.split()[0]]=''
                                    hap2_re_al[line.split()[0]]+=line.split()[1]

                        else:
                            continue
        
     
                        
                        seq_cons_nosub=''
                        seq_cons=''
                        seq_cons_nodel_2=''

                        seq_len=len(list(hap2_re_al.values())[0])

                        for i in range(seq_len):
                            counts={'A':0,'G':0,'T':0,'C':0,'-':0}
                            for read in hap2_re_al.values():
                                counts[read[i]]+=1
                            top_list=sorted(list(counts.items()), key=lambda x: x[1])
                            top=top_list[-1][0]
                            seq_cons_nosub+=top
                            if top=='-':
                                    if top_list[-2][1]/len(hap2_re_al.values())>=0.2:
                                        #print(top_list[-1][1],top_list[-2][1])
                                        top=top_list[-2][0]
                                        #print('subs %s for -' %top)
                                        #print(top_list)

                            seq_cons+=top
                            top ='' if top=='-' else top
                            seq_cons_nodel_2+=top


                        hap2_cons=seq_cons_nodel_2

                        ref=''.join(ref_df.reindex(range(v_pos-window_before,v_pos+window_after+1)).ref.to_list())

                        pair_1=pairwise(hap1_cons,ref)
                        pair_2=pairwise(hap2_cons,ref)
                        seq_ts_all,seq_ts=pair_1[0][0],pair_1[0][1]
                    
                        cnt=0
                        cnt_del=0

                        while cnt<window_before:
                            if seq_ts[cnt_del]!='-':
                                cnt+=1
                            cnt_del+=1

                        cnt=0
                        var_1=False
                        var_type_1=None
                        while cnt<10:
                            if cnt==6 and var_1==False:
                                break
                            if seq_ts[cnt_del+1+cnt]!=seq_ts_all[cnt_del+1+cnt] and (seq_ts[cnt_del+1+cnt]=='-' or seq_ts_all[cnt_del+1+cnt]=='-'):
                                var_1=True
                                if seq_ts[cnt_del+1+cnt]=='-':
                                    if var_type_1!='del':
                                        var_type_1='ins'
                                    else:
                                        var_1=False
                                        break
                                else:
                                    if var_type_1!='ins':
                                        var_type_1='del'
                                    else:
                                        var_1=False
                                        break
                            elif var_1==True:
                                break
                            cnt+=1

                        ref_out=''
                        alt_out=''
                        var_1_alleles=[]
                        if var_1:
                            ref_out=seq_ts[cnt_del:cnt_del+1+cnt].replace('-','')
                            alt_out=seq_ts_all[cnt_del:cnt_del+1+cnt].replace('-','')
                            if len(alt_out)>=1:
                                var_1_alleles=(v_pos,alt_out)

                        ref_pair1=ref_out
                        seq_ts_all,seq_ts=pair_2[0][0],pair_2[0][1]
                  
                        cnt=0
                        cnt_del=0
                        while cnt<window_before:
                            if seq_ts[cnt_del]!='-':
                                cnt+=1
                            cnt_del+=1

                        cnt=0
                        var_2=False
                        var_type_2=None

                        while cnt<10:
                            if cnt==6 and var_2==False:
                                break
                            if seq_ts[cnt_del+1+cnt]!=seq_ts_all[cnt_del+1+cnt] and (seq_ts[cnt_del+1+cnt]=='-' or seq_ts_all[cnt_del+1+cnt]=='-'):
                                    var_2=True
                                    if seq_ts[cnt_del+1+cnt]=='-':
                                        if var_type_2!='del':
                                            var_type_2='ins'
                                        else:
                                            var_2=False
                                            break
                                    else:
                                        if var_type_2!='ins':
                                            var_type_2='del'
                                        else:
                                            var_2=False
                                            break
                            elif var_2==True:
                                    break

                            cnt+=1

                        ref_out=''
                        alt_out=''

                        var_2_alleles=[]
                        if var_2:
                            ref_out=seq_ts[cnt_del:cnt_del+1+cnt].replace('-','')
                            alt_out=seq_ts_all[cnt_del:cnt_del+1+cnt].replace('-','')
                            if len(alt_out)>=1:
                                var_2_alleles=(v_pos,alt_out)

                        ref_pair2=ref_out        

                        if var_1_alleles and var_2_alleles:
                            if var_1_alleles[1]==var_2_alleles[1]:
                                res=(v_pos,(1,1),ref_pair1, var_1_alleles[1])
                                if len(res[2])>0 and len(res[3])>0:
                                    if var_1_alleles[-3:]!=var_1_alleles[-1]*3 and ref_pair1[-3:]!=ref_pair1[-1]*3:
                                        output_list.append(res)

                            else:
                                ref_pair=ref_pair1 if len(ref_pair1)>len(ref_pair2) else ref_pair2
                                res=(v_pos,(1,2),ref_pair,var_1_alleles[1], var_2_alleles[1])
                            
                                if len(res[2])>0 and len(res[3])>0 and len(res[4])>0:
                                    if var_1_alleles[-3:]!=var_1_alleles[-1]*3 and var_2_alleles[-3:]!=var_2_alleles[-1]*3 and ref_pair[-3:]!=ref_pair[-1]*3:
                                        output_list.append(res)

                        elif var_1_alleles:
                            res=(v_pos,(1,0),ref_pair1,var_1_alleles[1])
                            
                            if len(res[2])>0 and len(res[3])>0:
                                if var_1_alleles[-3:]!=var_1_alleles[-1]*3 and ref_pair1[-3:]!=ref_pair1[-1]*3:
                                    output_list.append(res)

                        elif var_2_alleles:
                            res=(v_pos,(0,1),ref_pair2,var_2_alleles[1])
                            if len(res[2])>0 and len(res[3])>0:
                                if var_2_alleles[-3:]!=var_2_alleles[-1]*3 and ref_pair2[-3:]!=ref_pair2[-1]*3:
                                    output_list.append(res)

                for pread in pcol.pileups:
                    if pread.query_position:
                        prev_cand_reads[pread.alignment.query_name]=pread.query_position
        except KeyError:
            continue
            
    print('%d done' %params['start'])
    return output_list
def generate_calls(params):
    #482375
    file_name=os.path.join(params['out_path'],'%s.%d-%d.vcf' %(params['chrom'],params['start'],params['end']))
    
    calls=indel_calls(params)

    
    f= open(file_name,'w')

    f.write('##fileformat=VCFv4.2\n')
    f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    c='##contig=<ID=%s>\n' %chrom
    f.write('##contig=<ID=%s>\n' %chrom)


    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write('##FORMAT=<ID=GP,Number=1,Type=Integer,Description="Genotype Probability">\n')
    f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
    for cand in calls:
            if cand[1]==(1,2):
                s='%s\t%d\t.\t%s\t%s,%s\t%d\t%s\t.\tGT:GP\t%s/%s:%d\n' %(chrom,cand[0], cand[2], cand[3], cand[4], 0,'PASS',cand[1][0],cand[1][1],0)
                f.write(s)
            else:
                s='%s\t%d\t.\t%s\t%s\t%d\t%s\t.\tGT:GP\t%s/%s:%d\n' %(chrom,cand[0], cand[2], cand[3], 0,'PASS',cand[1][0],cand[1][1],0)
                f.write(s)
    f.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #-r chromosome region   -m mode   -bam bam file   -ref reference file   -vcf ground truth variants   -o output path
    parser.add_argument("-chrom", "--chrom", help="Chromosome region")
    parser.add_argument("-bam", "--bam", help="Bam file")
    parser.add_argument("-ref", "--ref", help="Size")
    parser.add_argument("-vcf", "--vcf", help="Ground truth variants")
    parser.add_argument("-o", "--output", help="Output path")
    parser.add_argument("-cpu", "--cpu", help="CPUs",type=int)
    parser.add_argument("-bed", "--bed", help="BED file")
    parser.add_argument("-mincov", "--mincov", help="min coverage",type=int)
    parser.add_argument("-start", "--start", help="start",type=int)
    parser.add_argument("-end", "--end", help="end",type=int)
    parser.add_argument("-tree", "--tree", help="Haplotype tree path")

    
    args = parser.parse_args()
    
    chrom=args.chrom
    
    if not args.end:
        end=chrom_length[chrom]
    else:
        end=args.end
    
    if not args.start:
        start=1
        
    in_dict={'chrom':chrom,'start':args.start,'end':args.end, 'sam_path':args.bam, 'fasta_path':args.ref,\
             'out_path':args.output, 'cpu':args.cpu, 'bed':args.bed,'mincov':args.mincov,'tree_path':args.tree}    
    
    t=time.time()
    generate_calls(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
