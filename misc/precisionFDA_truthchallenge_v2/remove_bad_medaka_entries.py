from pysam import VariantFile
import sys

bcfin=VariantFile(sys.argv[1])
print(str(bcfin.header),end='')
for i in range(1,23):
    for rec in bcfin.fetch('chr%d' %i):
        if rec.ref not in rec.alts:
            print(str(rec),end='')
