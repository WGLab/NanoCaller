import sys,pysam, time,os,re,copy,argparse,gzip,itertools
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
from intervaltree import Interval, IntervalTree
from matplotlib import pyplot as plt

bcf_in=VariantFile(sys.argv[1])
with open('bad_nbr_sites_25.chr'+sys.argv[2],'w') as file:
    for c in [sys.argv[2]]:# range(1,23):
        nc_calls=[]
        chrom='chr%d' %int(c)
        for rec in bcf_in.fetch(chrom):
                nc_calls.append(rec.pos)
        nc_calls_df=pd.DataFrame(nc_calls)
        nc_calls_df.rename(columns={0:'pos'},inplace=True)
        nc_calls_df.set_index('pos',inplace=True,drop=False)

        nc_calls_tree=IntervalTree(Interval(pos-1000, pos+1000, "%d-%d" % (pos-1000, pos+1000)) for pos in nc_calls_df.pos)
        nc_calls_df[['nbr-nc_calls']]=nc_calls_df[['pos']].applymap(lambda x:len(nc_calls_tree[x]))

        new=nc_calls_df[nc_calls_df['nbr-nc_calls']>25].pos
        for x in new:
            file.write('%s\t%d\n' %(chrom,x))
        print('%s done' %chrom)
