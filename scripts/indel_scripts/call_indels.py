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

window_before,window_after=30,50

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

def msa_consensus(seq_list,ref):
    dnd_tmp_file,dnd_tmp= mkstemp()
    fa_tmp_file=''
    seq_count=0
    for b in seq_list:
        fa_tmp_file+='>hap1seq%d\n'%seq_count
        fa_tmp_file+= '%s\n' %b
        seq_count+=1
        
    msa_process = Popen(['t_coffee', '-infile=stdin','-outfile','stdout','-out_lib=/dev/null', '-no_warning', '-quiet','-newtree=%s'%dnd_tmp,'-type=dna'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    hap_file=msa_process.communicate(input=fa_tmp_file.encode('utf-8'))[0]

    realignments={}
    hap_file=hap_file.decode('utf-8')
    if 'hap' in hap_file:
        for line in hap_file.split('\n'):
            if 'hap' in line:
                if line.split()[0] not in realignments.keys():
                    realignments[line.split()[0]]=''
                realignments[line.split()[0]]+=line.split()[1]

    else:
        print('msa_process failed', flush=True)
        return (False,None,None,None)

    counts_df=pd.DataFrame([list(z) for z in realignments.values()]).apply(pd.value_counts).reindex(['A','G','T','C','-']).fillna(0)/len(seq_list)
    counts_df.loc['-']=counts_df.loc['-']-0.2
    seq_cons_nodel=''.join(counts_df.idxmax().to_list())

    pairwise_alignments=pairwise(seq_cons_nodel,ref)


    seq_ts_all,seq_ts=pairwise_alignments[0][0],pairwise_alignments[0][1]
    cnt=0
    cnt_del=0

    while cnt<window_before:
        if seq_ts[cnt_del]!='-':
            cnt+=1
        cnt_del+=1

    cnt=0
    var=False
    var_type=None
    while cnt<10:
        if seq_ts[cnt_del+1+cnt]!=seq_ts_all[cnt_del+1+cnt] and (seq_ts[cnt_del+1+cnt]=='-' or seq_ts_all[cnt_del+1+cnt]=='-'):
            var=True
            if seq_ts[cnt_del+1+cnt]=='-':
                if var_type!='del':
                    var_type='ins'
                else:
                    var=False
                    break
            else:
                if var_type!='ins':
                    var_type='del'
                else:
                    var=False
                    break
        elif var==True:
            break
        cnt+=1

    ref_out=''
    alt_out=''
    
    
    if var:
        ref_out=seq_ts[cnt_del:cnt_del+1+cnt].replace('-','')
        alt_out=seq_ts_all[cnt_del:cnt_del+1+cnt].replace('-','')
        
    if len(alt_out)==0:
        return (False,None,None,None)
    else:
        return (var,var_type,alt_out,ref_out)
    

def indel_calls(params):
    fasta_path=params['fasta_path']
    fastafile=pysam.FastaFile(fasta_path)
    sam_path=params['sam_path']
    samfile = pysam.Samfile(sam_path, "rb")
    chrom=params['chrom']
    bed_path=params['bed']

    hap_df=pd.read_csv('%s.%s.bam.info' %(params['hap_path'],chrom), delim_whitespace=True)

    hap_reads_0=set(hap_df[(hap_df.haplotype==0)]['#readname'])
    hap_reads_1=set(hap_df[(hap_df.haplotype==1)]['#readname'])
    

    if bed_path:
        bed_file=os.path.join(bed_path,'%s.bed' %params['chrom'])
        
        file=open(bed_file)
        content=[x.rstrip('\n') for x in file]
        file.close()
        
        content=[x.split('\t')[1:] for x in content]
        content=[(int(x[0]),int(x[1])) for x in content]
        bed_tree=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)


    output_list=[]

    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,params['start']-300),params['end']+300+1), fastafile.fetch(params['chrom'], max(1,params['start']-300)-1,params['end']+300)) }

    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)


    for pcol in samfile.pileup(chrom, params['start']-1, params['end'], min_base_quality=0, flag_filter=0x4|0x100|0x200|0x400|0x800, truncate=True):
        v_pos=pcol.pos+1
        if v_pos>=params['end']:
            pass#return output_list
        if bed_path:
            if not bed_tree[pcol.pos+1]:
                continue
        try:
            n=pcol.get_num_aligned()
            if n>params['mincov']:
                nbr=''.join(ref_df.reindex(range(pcol.pos-4,pcol.pos+6)).ref.to_list())
                check=repeat_check(nbr)
                
                var_0,var_1=False,False
                
                if check==False:                    
                    seq=''.join(pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True))
                    
                    name=pcol.get_query_names()
                    seq_0=''.join([base for base,read_name in zip(seq,name) if read_name in hap_reads_0])
                    p_count_0=seq_0.count('+')/len(seq_0) if len(seq_0)>0 else 0
                    n_count_0=(seq_0.count('-')+seq_0.count('*'))/len(seq_0) if len(seq_0)>0 else 0
                    
                    seq_1=''.join([base for base,read_name in zip(seq,name) if read_name in hap_reads_1])
                    p_count_1=seq_1.count('+')/len(seq_1) if len(seq_1)>0 else 0
                    n_count_1=(seq_1.count('-')+seq_1.count('*'))/len(seq_1) if len(seq_1)>0 else 0
                    
                    d={'hap0':{},'hap1':{}}
                    
                    if p_count_0>=0.3 or n_count_0>=0.6:
                        for pread in pcol.pileups:
                            if pread.alignment.qname in hap_reads_0:
                                d['hap0'][pread.alignment.qname]=pread.alignment.query_sequence[pread.query_position_or_next-window_before:pread.query_position_or_next+window_after]
                                
                        if len(d['hap0'])>=params['mincov']:
                            ref=''.join(ref_df.reindex(range(v_pos-window_before,v_pos+window_after+1)).ref.to_list())
                            var_0,var_type_0,alt_out_0,ref_out_0=msa_consensus(d['hap0'].values(),ref)
                            
                        
                    if p_count_1>=0.3 or n_count_1>=0.6:
                        
                        for pread in pcol.pileups:
                            if pread.alignment.qname in hap_reads_1:
                                d['hap1'][pread.alignment.qname]=pread.alignment.query_sequence[pread.query_position_or_next-window_before:pread.query_position_or_next+window_after]
                                
                        if len(d['hap1'])>=params['mincov']:
                            ref=''.join(ref_df.reindex(range(v_pos-window_before,v_pos+window_after+1)).ref.to_list())
                            var_1,var_type_1,alt_out_1,ref_out_1=msa_consensus(d['hap1'].values(),ref)
                            
                    v_pos=int(v_pos)        
                    if var_0 and var_1:
                        if alt_out_0==alt_out_1:
                            res=(v_pos,'1/1',ref_out_0, alt_out_0)
                            if len(ref_out_0)>0:
                                if alt_out_0[-3:]!=alt_out_0[-1]*3 and ref_out_0[-3:]!=ref_out_0[-1]*3:
                                    output_list.append(ref_out_0)

                        else:
                            ref_pair=ref_out_0 if len(ref_out_0)>len(ref_out_1) else ref_out_1
                            res=(v_pos,'1|2',ref_pair,alt_out_0, alt_out_1)

                            if len(ref_pair)>0:
                                if (alt_out_0[-3:]!=alt_out_0[-1]*3 or alt_out_1[-3:]!=alt_out_1[-1]*3) and ref_pair[-3:]!=ref_pair[-1]*3:
                                    output_list.append(res)

                    elif var_0:
                        res=(v_pos,'1|0',ref_out_0,alt_out_0)

                        if len(ref_out_0)>0:
                            if alt_out_0[-3:]!=alt_out_0[-1]*3 and ref_out_0[-3:]!=ref_out_0[-1]*3:
                                output_list.append(res)

                    elif var_1:
                        res=(v_pos,'0|1',ref_out_1,alt_out_1)
                        if len(ref_out_1)>0:
                            if alt_out_1[-3:]!=alt_out_1[-1]*3 and ref_out_1[-3:]!=ref_out_1[-1]*3:
                                output_list.append(res)


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
            if cand[1]=='1|2':
                s='%s\t%d\t.\t%s\t%s,%s\t0\tPASS\t.\tGT:GP\t%s:30\n' %(chrom,cand[0], cand[2], cand[3], cand[4],cand[1])
                f.write(s)
                
            else:
                s='%s\t%d\t.\t%s\t%s\t0\tPASS\t.\tGT:GP\t%s:30\n' %(chrom,cand[0], cand[2], cand[3],cand[1])
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
    parser.add_argument("-hap", "--hap", help="Haplotype information file")

    
    args = parser.parse_args()
    
    chrom=args.chrom
    
    if not args.hap:
        hap='/home/ahsanm/test_nano/ts_HG2/tr_HG1/phasing/HG002'
    
    if not args.end:
        end=chrom_length[chrom]
    else:
        end=args.end
    
    if not args.start:
        start=1
    
    in_dict={'chrom':chrom,'start':args.start,'end':args.end, 'sam_path':args.bam, 'fasta_path':args.ref,\
             'out_path':args.output, 'cpu':args.cpu, 'bed':args.bed,'mincov':args.mincov,'hap_path':hap}    
    
    t=time.time()
    generate_calls(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
