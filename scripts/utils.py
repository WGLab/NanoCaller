import sys,gzip
import pandas as pd
import numpy as np

from matplotlib import pyplot as plt


def print_vcf(file,chrom,df):
    with open(file,'w') as f:
        
        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##contig=<ID=%s>\n' %chrom)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus Genotype across all datasets with called genotype">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
        for v in df.values:
            s='%s\t%d\t.\t%s\t%s\t%d\tPASS\t.\tGT\t%s\n' %(chrom,v[0],v[1],v[2],0,v[3])
            f.write(s)

def read_pileups_from_file(fname,dims,mode):
    lines={}
    if mode=='train':
        with open(fname,'rb') as file:
            for l in file:
                l=l.decode('utf-8')[:-1]
                pos,gtype,allele,ref,m1,m2=l.split(':')
                mm1=np.array(list(m1)).astype(np.int8).reshape((dims[0],dims[1],dims[2]-1))
                m2=m2.replace(".", "")
                mm2=np.array(m2.split('|')).astype(np.int8).reshape((dims[0],dims[1],1))
                p_mat=np.dstack((mm1,mm2))

                lines[pos]=(int(pos),p_mat,int(gtype[0]),int(allele),int(ref))
    
    else:
        with open(fname,'r') as file:
            for l in file:
                l=l[:-1]
                pos,ref,m1,m2=l.split(':')
                mm1=np.array(list(m1)).astype(np.int8).reshape((dims[0],dims[1],dims[2]-1))
                m2=m2.replace(".", "")
                mm2=np.array(m2.split('|')).astype(np.int8).reshape((dims[0],dims[1],1))
                p_mat=np.dstack((mm1,mm2))

                lines[pos]=(int(pos),p_mat,int(ref))
    
    return lines

def pileup_image_from_mat(in_mat):
        t=np.copy(in_mat[:,:,-1])
        t[t>0]=-1
        t[t==0]=1
        t[t==-1]=0

        new=np.dstack([in_mat[:,:,:5],t])
        p_mat=np.argmax(new,2)
        
        #deletion=purple,  1A=green,   2G=orange,    3T=red,    4C=blue
        color={5:np.array([[[0,0,0]]]),4:np.array([[[255, 255, 255]]]),0:np.array([[[0,   255, 0]]]),1:np.array([[[206, 132, 47]]]),2:np.array([[[255,   0, 0]]]),3:np.array([[[0,   0, 255]]])}
        data = np.zeros((p_mat.shape[0],p_mat.shape[1], 3)).astype(int)
        for j in range(6):
            new=(p_mat==j).astype(int)
            data+=np.repeat(new[:, :, np.newaxis], 3, axis=2)*color[j]
        
        plt.figure(figsize=(20,10))
        plt.imshow(data)
        plt.show()

def plot_training_stats(stats):
        plt.plot(range(stats.shape[0]), stats[:,0], 'b', label='Training loss')
        plt.plot(range(stats.shape[0]), stats[:,1], 'r', label='Test loss')
        plt.title('Training and Test loss')
        plt.xlabel('Epochs ',fontsize=16)
        plt.ylabel('Loss',fontsize=16)
        plt.legend()
        plt.figure()
        plt.show()
        fig1 = plt.gcf()
        fig1.savefig('loss.png',dpi=1000)

        plt.plot(range(stats.shape[0]), stats[:,2], 'b', label='Training Accuracy')
        plt.plot(range(stats.shape[0]), stats[:,3], 'r', label='Test Accuracy')
        plt.title('Training and Test Accuracy')
        plt.xlabel('Epochs ',fontsize=16)
        plt.ylabel('Accuracy',fontsize=16)
        plt.legend()
        plt.figure()
        plt.show()
        fig2 = plt.gcf()
        fig2.savefig('accuracy.png',dpi=1000)