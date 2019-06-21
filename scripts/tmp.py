import sys,pysam, time,os,re,copy,argparse,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
import generate_pileups as gcp
from matplotlib import pyplot as plt
from importlib import reload
import utils


fasta_path='/home/ahsanm1/shared_data/Ref_Genome/fasta/hg19/hg19.fa'
vcf_path='/home/ahsanm1/shared_data/NA12878/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz '
amplicon_path='/home/ahsanm1/umair_wlab/data/CYP2D6/CYP2D6.bam'

t=time.time()
amp_params={'mode':'testing','threshold':0.25,'bed':False,'chrom':'chr22','start':42522077,'end':42527144,'sam_path':amplicon_path,'out_path':'/home/ahsanm1/umair_wlab/data/CYP2D6/1','fasta_path':fasta_path,'vcf_path':vcf_path,'depth':32,'window':16,'cpu':10}
gcp.generate(amp_params)

elapsed=time.time()-t
print('Elapsed Time %.2f' %elapsed)

t=time.time()
amp_params={'mode':'testing','threshold':0.25,'bed':False,'chrom':'chr22','start':42522077,'end':42527144,'sam_path':amplicon_path,'out_path':'/home/ahsanm1/umair_wlab/data/CYP2D6/2','fasta_path':fasta_path,'vcf_path':vcf_path,'depth':64,'window':16,'cpu':10}
gcp.generate(amp_params)

elapsed=time.time()-t
print('Elapsed Time %.2f' %elapsed)

t=time.time()
amp_params={'mode':'testing','threshold':0.25,'bed':False,'chrom':'chr22','start':42522077,'end':42524144,'sam_path':amplicon_path,'out_path':'/home/ahsanm1/umair_wlab/data/CYP2D6/3','fasta_path':fasta_path,'vcf_path':vcf_path,'depth':128,'window':16,'cpu':10}
gcp.generate(amp_params)

elapsed=time.time()-t
print('Elapsed Time %.2f' %elapsed)

t=time.time()
amp_params={'mode':'testing','threshold':0.25,'bed':False,'chrom':'chr22','start':42522077,'end':42527144,'sam_path':amplicon_path,'out_path':'/home/ahsanm1/umair_wlab/data/CYP2D6/4','fasta_path':fasta_path,'vcf_path':vcf_path,'depth':256,'window':16,'cpu':10}
gcp.generate(amp_params)

elapsed=time.time()-t
print('Elapsed Time %.2f' %elapsed)

t=time.time()
amp_params={'mode':'testing','threshold':0.25,'bed':False,'chrom':'chr22','start':42522077,'end':42527144,'sam_path':amplicon_path,'out_path':'/home/ahsanm1/umair_wlab/data/CYP2D6/5','fasta_path':fasta_path,'vcf_path':vcf_path,'depth':512,'window':16,'cpu':10}
gcp.generate(amp_params)

elapsed=time.time()-t
print('Elapsed Time %.2f' %elapsed)

t=time.time()
amp_params={'mode':'testing','threshold':0.25,'bed':False,'chrom':'chr22','start':42522077,'end':42527144,'sam_path':amplicon_path,'out_path':'/home/ahsanm1/umair_wlab/data/CYP2D6/6','fasta_path':fasta_path,'vcf_path':vcf_path,'depth':1024,'window':16,'cpu':10}
gcp.generate(amp_params)

elapsed=time.time()-t
print('Elapsed Time %.2f' %elapsed)

