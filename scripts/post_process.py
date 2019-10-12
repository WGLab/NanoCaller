import sys,pysam, time,os,re,copy,argparse,gzip,itertools
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree


chrom_length={'chr1':248956422, 'chr2':242193529, 'chr3':198295559, 'chr4':190214555, 'chr5':181538259, 'chr6':170805979, \
             'chr7':159345973, 'chr8':145138636, 'chr9':138394717, 'chr10':133797422, 'chr11':135086622, 'chr12':133275309,\
             'chr13':114364328, 'chr14':107043718, 'chr15':101991189, 'chr16':90338345, 'chr17':83257441, 'chr18':80373285,\
             'chr19':58617616, 'chr20':64444167, 'chr21':46709983, 'chr22':50818468, 'chrX':156040895, 'chrY':57227415}

mapping={'A':0,'G':1,'T':2,'C':3,'*':4,'N':4}
    
def create_pileup_df(params):
    
    with open(os.path.join(params['tree_path'],'%s.tree' %params['chrom']),'r') as file:
        content=file.readlines()
    
    hap_tree=IntervalTree()

    for line in content:
        line=line.rstrip('\n')
        line=line.split('>')
        begin,end=int(line[0].split('-')[0]),int(line[0].split('-')[1])
        hap_tree[begin:end+1]=set(line[1].split(':'))
        
        
    bcf_in = VariantFile(params['vcf_path'])
    neg_vcf=params['vcf_path'][:-3]+'.neg'
    fastafile=pysam.FastaFile(params['fasta_path'])
    samfile = pysam.Samfile(params['sam_path'], "rb")
    bed_path=params['bed']
    
    neg_cand=pd.read_csv(neg_vcf,header=None)
    neg_cand=neg_cand[(neg_cand[0]>=params['start'])&(neg_cand[0]<=params['end'])]
    neg_cand=set(neg_cand[0])

    
    cand=[]
    
    start=max(params['start']-100000,0)
    end=min(params['end']+100000,chrom_length[params['chrom']])
    
    tot_cand=[]
    for rec in bcf_in.fetch(params['chrom'],start=start,end=end):
        tot_cand.append(rec.pos)

    tot_cand=set(tot_cand)
    
    mapping={'A':0,'G':1,'T':2,'C':3,'*':4,'N':4}
    
    pileup_dict={}
    
    if bed_path:
        bed_file=os.path.join(bed_path,'%s.bed' %params['chrom'])
        with open(bed_file) as file:
            content=[x.rstrip('\n') for x in file]

        content=[x.split('\t')[1:] for x in content]
        content=[(int(x[0]),int(x[1])) for x in content]
        bed_tree=IntervalTree(Interval(begin, end, "%d-%d" % (begin, end)) for begin, end in content)
        
    ref_dict={j:s.upper() if s in 'AGTC' else '*' for j,s in zip(range(max(1,start-40),end+40+1), fastafile.fetch(params['chrom'], max(1,start-40)-1,end+40)) }
    
    ref_df=pd.DataFrame(list(ref_dict.items()), columns=['pos', 'ref'])
    ref_df.set_index('pos',drop=False,inplace=True)
        
    output=[]
    
    for pcol in samfile.pileup(chrom,max(0,start-1-30),end+30,min_base_quality=0,\
                                           flag_filter=0x4|0x100|0x200|0x400|0x800,truncate=True):
            v_pos=pcol.pos+1
            cand_ref=ref_df.loc[v_pos].ref
            
            if bed_path:
                if not bed_tree[v_pos]:
                    continue
            
            if cand_ref!='*' and hap_tree[v_pos] and v_pos not in tot_cand:
                seq=''.join([x[0] for x in pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                #seq_nodel=''.join([pcol.get_query_sequences( mark_matches=False, mark_ends=False,add_indels=True)]).upper()
                n=pcol.get_num_aligned()
                #dp_nodel=len(seq_nodel)
                


                if n>=params['mincov'] and v_pos>=start and v_pos<=end:
                    alt_freq=max([x[1] for x in Counter(seq).items() if (x[0]!=cand_ref and x[0] in 'AGTC')]+[0])/n
                    
                    if 0.05<=alt_freq or v_pos in neg_cand:
                        cand_info={n:s for n,s in zip(pcol.get_query_names(),seq)}
                        common_reads=[]
                        for iv in hap_tree[v_pos]:
                            common_reads.append(set(iv.data)&set(cand_info.keys()))

                        common_reads=sorted(common_reads, key=lambda x: len(x),reverse=True)

                        #print(cnd_df.loc[v_pos])
                        #print('\n')
                        check=[]
                                                
                        if len(common_reads)<2:
                            continue
                            
                        hap1=cand_ref
                        hap2=cand_ref
                        rd_gp =common_reads[0]
                        
                        seq1=''.join([cand_info[x] for x in rd_gp])
                        ref_count1=seq1.count(cand_ref)
                        cc1=max([(s,seq1.count(s)) for s in 'AGTC' if s!=cand_ref], key=lambda x:x[1])                           
                        
                        
                        '''if v_pos in neg_cand:
                            if (cc1[1]/len(seq1)>=0.4 or cc1[1]>=ref_count1-2) and  cc1[1]>=2:
                                hap1=cc1[0]'''

                        
                        if cc1[1]>=ref_count1-1  and  cc1[1]>=4:
                            hap1=cc1[0]

                        
                        rd_gp =common_reads[1]
                        
                        seq2=''.join([cand_info[x] for x in rd_gp])
                        ref_count2=seq2.count(cand_ref)
                        cc2=max([(s,seq2.count(s)) for s in 'AGTC' if s!=cand_ref], key=lambda x:x[1])
                            
                        '''if v_pos in neg_cand:
                            if (cc2[1]/len(seq2)>=0.4 or cc2[1]>=ref_count2-2) and  cc2[1]>=2:
                                hap2=cc2[0]'''

                        
                        if cc2[1]>=ref_count2-1 and  cc2[1]>=4:
                            hap2=cc2[0]

                            
                        ts_tv={'A':2,'G':3,'T':5,'C':7}
                        
                        if hap1==cand_ref and hap2!=cand_ref:
                            a=cc2[1]/len(seq2)
                            output.append((v_pos,(0,1),cand_ref,hap2,int(min(999,-100*np.log10(1e-10+1-a)))))

                        elif hap2==cand_ref and hap1!=cand_ref:
                            a=cc1[1]/len(seq1)
                            output.append((v_pos,(0,1),cand_ref,hap1,int(min(999,-100*np.log10(1e-10+1-a)))))
                       

                        elif hap1!=cand_ref and hap2!=cand_ref and hap1==hap2:
                                a=(cc1[1]+cc2[1])/(len(seq2)+len(seq1))
                                output.append((v_pos,(1,1),cand_ref,hap1,int(min(999,-100*np.log10(1e-10+1-a)))))
                                
                        elif hap1!=cand_ref and hap2!=cand_ref and hap1!=hap2:
                                a=(cc1[1]+cc2[1])/(len(seq2)+len(seq1))
                                output.append((v_pos,(1,2),cand_ref,hap1,hap2,int(min(999,-100*np.log10(1e-10+1-a)))))
 
    print('%d-%d calling done'%(params['start'],params['end']),flush=True)
    return output

def generate_calls(params):
    bcf_in = VariantFile(params['vcf_path'])
    fastafile=pysam.FastaFile(params['fasta_path'])
    samfile = pysam.Samfile(params['sam_path'], "rb")

    dict_list=[]

    for i in range(0,chrom_length[params['chrom']],int(1e6)):
        tmp_dict=copy.deepcopy(params)
        tmp_dict['start']=i
        tmp_dict['end']=i+1000000
        dict_list.append(tmp_dict)

    pool = mp.Pool(processes=params['cpu'])
    calls = pool.map(create_pileup_df, dict_list)

    
    with open(params['out_path'],'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        

        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GP,Number=1,Type=Integer,Description="Genotype Probability">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
        for output in calls:
            for cand in output:
 
                if cand[1]==(1,2):
                    s='%s\t%d\t.\t%s\t%s,%s\t%d\t%s\t.\tGT:GP\t%s/%s:%d\n' %(chrom,cand[0], cand[2], cand[3], cand[4], 0,'PASS',cand[1][0],cand[1][1],cand[-1])
                    f.write(s)
                else:
                    s='%s\t%d\t.\t%s\t%s\t%d\t%s\t.\tGT:GP\t%s/%s:%d\n' %(chrom,cand[0], cand[2], cand[3], 0,'PASS',cand[1][0],cand[1][1],cand[-1])
                    f.write(s)

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
        
    in_dict={'chrom':chrom,'start':start,'end':end, 'sam_path':args.bam, 'fasta_path':args.ref, 'vcf_path':args.vcf,\
             'out_path':args.output, 'cpu':args.cpu, 'bed':args.bed,'mincov':args.mincov,'tree_path':args.tree}    
    
    t=time.time()
    generate_calls(in_dict)
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)