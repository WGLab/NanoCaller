from pysam import VariantFile
import sys

bcfin=VariantFile(sys.argv[1])
var_calls=[]
min_score=10000
max_score=-1
for i in range(1,23):
    for rec in bcfin.fetch('chr%d' %i):
        var_calls.append(rec)
        
        if rec.qual< min_score:
            min_score=rec.qual
        if rec.qual> max_score:
            max_score=rec.qual
            
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
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s''' %sys.argv[2]



print(header)
for rec in var_calls:
    gq=100*(rec.qual-min_score)/(max_score-min_score)+1
    rec.qual=gq
    s=str(rec).split('\t')
    print('%s\tGT:GQ\t%s:%.2f' %('\t'.join(s[:-2]),s[-1].split(':')[0],gq))