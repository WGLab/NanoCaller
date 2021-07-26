from pysam import VariantFile
import sys,random

bcfin=VariantFile(sys.argv[1])
var_calls=[]
for i in range(1,23):
    for rec in bcfin.fetch('chr%d' %i):
        var_calls.append(rec)
        
        
header='''##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##INFO=<ID=votes,Number=1,Type=Integer,Description="votes">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s''' %sys.argv[2]

print(header)
i=0
for original_rec in var_calls:
    rec=original_rec.copy()
    base_counter={('N','N'):0,('A', 'A'): 0, ('A', 'C'): 0, ('A', 'G'): 0, ('A', 'T'): 0, ('C', 'C'): 0, ('C', 'G'): 0, ('C', 'T'): 0, ('G', 'G'): 0, ('G', 'T'): 0, ('T', 'T'): 0}
    score_counter={('N','N'):0,('A', 'A'): 0, ('A', 'C'): 0, ('A', 'G'): 0, ('A', 'T'): 0, ('C', 'C'): 0, ('C', 'G'): 0, ('C', 'T'): 0, ('G', 'G'): 0, ('G', 'T'): 0, ('T', 'T'): 0}
    for x in rec.samples.keys():
        gt=rec.samples[x]['GT']
        if gt==(None,None) or gt==(0,0):
            base_counter[('N','N')]+=1
        else:
            base_counter[tuple (sorted( [rec.alleles[gt[0]],rec.alleles[gt[1]]]))]+=1
            score_counter[tuple (sorted( [rec.alleles[gt[0]],rec.alleles[gt[1]]]))]+=rec.qual

    max_keys=[x for x in base_counter.keys() if base_counter[x]==max(base_counter.values())]

    if sys.argv[3]=='v1':
        if (len(max_keys)==1 and ('N','N') in max_keys) or max(base_counter.values())<=len(rec.samples.keys())//2:
            continue
    else:
        if (len(max_keys)==1 and ('N','N') in max_keys):
            continue
    max_keys=[x for x in max_keys if x!=('N','N')]
    output=random.choice(max_keys)
    if output[0]==output[1]:
        rec.alts=[output[0]]
        gt='1/1'

    else:
        if output[0]==rec.ref:
            rec.alts=[output[1]]
            gt='0/1'
        elif output[1]==rec.ref:
            rec.alts=[output[0]]
            gt='0/1'
        else:
            rec.alts=output
            gt='1/2'

    rec.qual=score_counter[output]/len(rec.samples.keys())
    s=str(rec).rstrip('\n').split('\t')
    print('\t'.join(s[:7])+'\tvotes=%d\tGT:GQ\t%s:%.2f' %(base_counter[output],gt,rec.qual))
