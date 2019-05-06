import sys,gzip
import pandas as pd
import numpy as np

from matplotlib import pyplot as plt




def read_pileups_from_file(fname):
    lines={}
    with gzip.open(fname,'rb') as file:
        for l in file:
            l=l.decode('utf-8')[:-1]
            [pos,gtype,rnames,m1,m2]=l.split(':')
            rnames=rnames.split('|')
            mm1=np.array(list(m1)).astype(np.int8).reshape((32,101,7))
            m2=m2.replace(".", "")
            mm2=np.array(m2.split('|')).astype(np.int8).reshape((32,101,1))
            p_mat=np.dstack((mm1,mm2))

            lines[pos]=(int(pos),int(gtype[0]),rnames,p_mat)
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