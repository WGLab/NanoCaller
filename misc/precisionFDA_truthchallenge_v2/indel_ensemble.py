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
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s''' %sys.argv[2]
gt_map={(0, 0): 0, (0, 1): 1, (0, 2): 1, (0, 3): 1, (0, 4): 1, (0, 5): 1, (0, 6): 1, (0, 7): 1, (0, 8): 1, (0, 9): 1, (1, 0): 1, (1, 1): 0, (1, 2): 1, (1, 3): 1, (1, 4): 1, (1, 5): 1, (1, 6): 1, (1, 7): 1, (1, 8): 1, (1, 9): 1, (2, 0): 1, (2, 1): 1, (2, 2): 0, (2, 3): 1, (2, 4): 1, (2, 5): 1, (2, 6): 1, (2, 7): 1, (2, 8): 1, (2, 9): 1, (3, 0): 1, (3, 1): 1, (3, 2): 1, (3, 3): 0, (3, 4): 1, (3, 5): 1, (3, 6): 1, (3, 7): 1, (3, 8): 1, (3, 9): 1, (4, 0): 1, (4, 1): 1, (4, 2): 1, (4, 3): 1, (4, 4): 0, (4, 5): 1, (4, 6): 1, (4, 7): 1, (4, 8): 1, (4, 9): 1, (5, 0): 1, (5, 1): 1, (5, 2): 1, (5, 3): 1, (5, 4): 1, (5, 5): 0, (5, 6): 1, (5, 7): 1, (5, 8): 1, (5, 9): 1, (6, 0): 1, (6, 1): 1, (6, 2): 1, (6, 3): 1, (6, 4): 1, (6, 5): 1, (6, 6): 0, (6, 7): 1, (6, 8): 1, (6, 9): 1, (7, 0): 1, (7, 1): 1, (7, 2): 1, (7, 3): 1, (7, 4): 1, (7, 5): 1, (7, 6): 1, (7, 7): 0, (7, 8): 1, (7, 9): 1, (8, 0): 1, (8, 1): 1, (8, 2): 1, (8, 3): 1, (8, 4): 1, (8, 5): 1, (8, 6): 1, (8, 7): 1, (8, 8): 0, (8, 9): 1, (9, 0): 1, (9, 1): 1, (9, 2): 1, (9, 3): 1, (9, 4): 1, (9, 5): 1, (9, 6): 1, (9, 7): 1, (9, 8): 1, (9, 9): 0}
print(header)
i=0
prev_len=0
prev_pos=0
prev_chrom='none'
for original_rec in var_calls:
    if original_rec.contig!=prev_chrom:
        prev_len=0
        prev_pos=0
    if original_rec.pos<=prev_pos+prev_len:
        continue
    rec=original_rec.copy()
    
    prev_chrom=original_rec.contig
    gt=[tuple(sorted(rec.samples[x]['GT'])) for x in rec.samples.keys() if rec.samples[x]['GT']!=(None,None)]
    caller=[x for x in rec.samples.keys() if rec.samples[x]['GT']!=(None,None)]
    if rec.samples['medaka']['GT']!=(None,None):
        gq=0
        for x in rec.samples.keys():
            if type(rec.samples[x]['GQ'])==tuple:
                if rec.samples[x]['GQ'][0]!=None:
                    gq+=rec.samples[x]['GQ'][0]
            else:
                if rec.samples[x]['GQ']!=None:
                    gq+=rec.samples[x]['GQ']
        tmp_gt=rec.samples['medaka']['GT']
        gt_type=gt_map[tmp_gt]
        rec.alts=set([rec.alleles[x] for x in rec.samples['medaka']['GT'] if rec.alleles[x]!=rec.ref])
        
        if len(rec.alts)==1:
            if gt_type==0:
                tmp_gt=(1,1)
            else:
                tmp_gt=(0,1)
        else:
            tmp_gt=(1,2)

        s=str(rec).rstrip('\n').split('\t')
        print('\t'.join(s[:-6])+'\tPASS\t.\tGT\t%d/%d' %(tmp_gt[0],tmp_gt[1]))
        prev_len=max([len(x) for x in original_rec.alleles])
        prev_pos=rec.pos
        
    elif rec.samples['clair']['GT']!=(None,None):
        gq=0
        for x in rec.samples.keys():
            if type(rec.samples[x]['GQ'])==tuple:
                if rec.samples[x]['GQ'][0]!=None:
                    gq+=rec.samples[x]['GQ'][0]
            else:
                if rec.samples[x]['GQ']!=None:
                    gq+=rec.samples[x]['GQ']
                    
        tmp_gt=rec.samples['clair']['GT']
        gt_type=gt_map[tmp_gt]
        rec.alts=set([rec.alleles[x] for x in rec.samples['clair']['GT'] if rec.alleles[x]!=rec.ref])
        
        if len(rec.alts)==1:
            if gt_type==0:
                tmp_gt=(1,1)
            else:
                tmp_gt=(0,1)
        else:
            tmp_gt=(1,2)

        s=str(rec).rstrip('\n').split('\t')
        print('\t'.join(s[:-6])+'\tPASS\t.\tGT\t%d/%d' %(tmp_gt[0],tmp_gt[1]))
        prev_len=max([len(x) for x in original_rec.alleles])
        prev_pos=rec.pos        
        
    elif rec.samples['nanocaller']['GT']!=(None,None):
        gq=0
        for x in rec.samples.keys():
            if type(rec.samples[x]['GQ'])==tuple:
                if rec.samples[x]['GQ'][0]!=None:
                    gq+=rec.samples[x]['GQ'][0]
            else:
                if rec.samples[x]['GQ']!=None:
                    gq+=rec.samples[x]['GQ']
                    
        tmp_gt=rec.samples['nanocaller']['GT']
        gt_type=gt_map[tmp_gt]
        rec.alts=set([rec.alleles[x] for x in rec.samples['nanocaller']['GT'] if rec.alleles[x]!=rec.ref])
        
        if len(rec.alts)==1:
            if gt_type==0:
                tmp_gt=(1,1)
            else:
                tmp_gt=(0,1)
        else:
            tmp_gt=(1,2)

        s=str(rec).rstrip('\n').split('\t')
        print('\t'.join(s[:-6])+'\tPASS\t.\tGT\t%d/%d' %(tmp_gt[0],tmp_gt[1]))
        prev_len=max([len(x) for x in original_rec.alleles])
        prev_pos=rec.pos       