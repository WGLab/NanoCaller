import sys,pysam, time,os,re,copy,argparse, tensorboard,gzip
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
import tensorflow as tf
import generate_pileups as gcp
from matplotlib import pyplot as plt
from importlib import reload
import CNN_model as model
import utils

test_path='/home/ahsanm1/umair_wlab/data/CYP2D6/1/chr22_pileups_test'
fasta_path='/home/ahsanm1/shared_data/Ref_Genome/fasta/hg19/hg19.fa'
chrom='chr22'
vcf_path='/home/ahsanm1/scratch/amp.vcf'
dims=[32,33,5]
cpu=4

params={'chrom':chrom,'cpu':cpu,'vcf_path':vcf_path,'test_path':test_path,'dims':dims,'model':'/home/ahsanm1/umair_wlab/data/NanoVar_data/models/june4/model-11'}

model.test_model(params)
