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


f_path='/home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/training/chr20/chr20_pileups_'
cpu=16
_,x_train,y_train,train_allele,train_ref= model.get_data(f_path+'pos',cpu=cpu)

negative_variants=[model.get_data(f_path+'neg_%d' %freq,cpu=cpu) for freq in [0,5,10,15,20,25]]

nx_train=np.vstack([tmp[1] for tmp in negative_variants])
ny_train=np.vstack([tmp[2] for tmp in negative_variants])
ntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
ntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

perm=np.random.permutation(len(nx_train))

np.take(nx_train,perm,axis=0,out=nx_train)
np.take(ny_train,perm,axis=0,out=ny_train)
np.take(ntrain_allele,perm,axis=0,out=ntrain_allele)
np.take(ntrain_ref,perm,axis=0,out=ntrain_ref)

f_path='/home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/training/chr19/chr19_pileups_'
_,vpx_train,vpy_train,vptrain_allele,vptrain_ref= model.get_data(f_path+'pos',cpu=cpu)

negative_variants=[model.get_data(f_path+'neg_%d' %freq,cpu=cpu) for freq in [0,5,10,15,20,25]]

vnx_train=np.vstack([tmp[1] for tmp in negative_variants])
vny_train=np.vstack([tmp[2] for tmp in negative_variants])
vntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
vntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

perm=np.random.permutation(len(vnx_train))

np.take(vnx_train,perm,axis=0,out=vnx_train)
np.take(vny_train,perm,axis=0,out=vny_train)
np.take(vntrain_allele,perm,axis=0,out=vntrain_allele)
np.take(vntrain_ref,perm,axis=0,out=vntrain_ref)

test_size=len(vnx_train)
vx_test,vy_test,vtest_allele,vtest_ref=np.vstack([vpx_train,vnx_train[:test_size]]),\
np.vstack([vpy_train,vny_train[:test_size]]),np.vstack([vptrain_allele,vntrain_allele[:test_size]]),\
np.vstack([vptrain_ref,vntrain_ref[:test_size]])

params={'cpu':24,'rate':0.001,'iters':100,'size':50,'dims': [32,33,5],'model':'/home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/training/model'}
data={'train_data':[(x_train,y_train,train_allele,train_ref),(nx_train,ny_train,ntrain_allele,ntrain_ref)],'test_data':(vx_test,vy_test,vtest_allele,vtest_ref)}

res=model.genotype_caller(params,input_type=data,data=data)