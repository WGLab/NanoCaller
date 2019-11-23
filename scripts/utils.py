import sys,os,psutil,subprocess,time
import pandas as pd
import numpy as np
import multiprocessing as mp

from matplotlib import pyplot as plt

def f1(x,y,z):
    a=z/(z+y)
    b=z/(z+x)
    print('Precision: %.4f    Recall: %.4f    F1: %.4f' %(a,b,2*a*b/(a+b)))
    
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

def get_train_test(params,mode='train',verbose=False,tot_list=None,i=None,path=None,skip_pos=False):
    cpu=params['cpu']
    n_input=params['dims']
    if mode=='train':

        if path:
            f_path=path
        else:
            f_path=params['train_path']
        
        _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos',cpu=cpu,dims=n_input)
        
        _,nx_train,ny_train,ntrain_allele,ntrain_ref=get_data(f_path+'neg.combined',cpu=cpu,dims=n_input)

        return (_,x_train,y_train,train_allele,train_ref), (_,nx_train,ny_train,ntrain_allele,ntrain_ref)

    else:
        cpu=params['cpu']
        n_input=params['dims']
        
        if path:
            f_path=path
        else:
            f_path=params['test_path']
        
        _,vpx_train,vpy_train,vptrain_allele,vptrain_ref= get_data(f_path+'pos',cpu=cpu,dims=n_input)
        _,vnx_test,vny_test,vntest_allele,vntest_ref=get_data(f_path+'neg.combined',cpu=cpu,dims=n_input)
        vx_test,vy_test,vtest_allele,vtest_ref =np.vstack([vpx_train,vnx_test]), np.vstack([vpy_train,vny_test]), np.vstack([vptrain_allele,vntest_allele]), np.vstack([vptrain_ref,vntest_ref])

        
        return (vx_test,vy_test,vtest_allele,vtest_ref)
        
            
def get_data(fname,a=None, b=None,dims=(32,33,5), cpu=4,mode='train',verbose=False):
    t=time.time()
    l=os.stat(fname).st_size
    
    if mode=='train':
        rec_size=14+dims[0]*dims[1]*dims[2]*6
        if a!=None and b!=None:
            my_array=[(fname,x,mode,dims) for x in range(a,b,1000*rec_size)]
        else:
            my_array=[(fname,x,mode,dims) for x in range(0,l,1000*rec_size)]
    else:
        rec_size=12+8+dims[0]*dims[1]*dims[2]*6
        if a!=None and b!=None:
            my_array=[(fname,x,mode,dims) for x in range(a,b,1000*rec_size)]
        else:
            my_array=[(fname,x,mode,dims) for x in range(0,l,1000*rec_size)]

    cpu=min(cpu,len(my_array))
    pool = mp.Pool(processes=cpu)
    results = pool.map(read_pileups_from_file, my_array)
    pool.close()  
    pool.join() 
    
    
    pos=np.vstack([res[0][:,np.newaxis] for res in results])
    mat=np.vstack([res[1] for res in results])
    
    '''if window:
            w=window
            l=dims[1]//2
            mat=mat[:,:,l-w:l+w+1,:]'''
            
    ref=np.vstack([res[2] for res in results])
    allele,gt=None,None
    if mode=='train':
        allele=np.vstack([res[3] for res in results])
        gt=np.vstack([res[4] for res in results])
    
        elapsed=time.time()-t

        if verbose:
            print('I/O Time Elapsed: %.2f seconds' %elapsed, flush = True)


        return pos,mat,gt,allele,ref
    
    else:
        dp=np.vstack([res[3][:,np.newaxis] for res in results])

        freq=np.vstack([res[4][:,np.newaxis] for res in results])
    
        elapsed=time.time()-t

        if verbose:
            print('I/O Time Elapsed: %.2f seconds' %elapsed, flush = True)

        return pos,mat,ref,dp,freq

def get_train_test_old(params,mode='train',verbose=False,tot_list=None,i=None,path=None,skip_pos=False):
    cpu=params['cpu']
    n_input=params['dims']
    if mode=='train':

        if path:
            f_path=path
        else:
            f_path=params['train_path']
        
        _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos',verbose=verbose,cpu=cpu,dims=n_input)
        
        negative_variants=[get_data(f_path+'neg.%d' %freq,cpu=cpu,verbose=verbose,dims=n_input) for freq in [0,5,10,15,20,25]]

        nx_train=np.vstack([tmp[1] for tmp in negative_variants])
        ny_train=np.vstack([tmp[2] for tmp in negative_variants])
        ntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
        ntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

        perm=np.random.permutation(len(nx_train))

        np.take(nx_train,perm,axis=0,out=nx_train)
        np.take(ny_train,perm,axis=0,out=ny_train)
        np.take(ntrain_allele,perm,axis=0,out=ntrain_allele)
        np.take(ntrain_ref,perm,axis=0,out=ntrain_ref)

        return (_,x_train,y_train,train_allele,train_ref), (_,nx_train,ny_train,ntrain_allele,ntrain_ref)
    
    elif mode=='train_chunk':
        cpu=params['cpu']
        n_input=params['dims']
        
        if path:
            f_path=path
        else:
            f_path=params['train_path']
        if skip_pos:
            _,x_train,y_train,train_allele,train_ref=None,None,None,None,None
        
        else:
            _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos',a=tot_list['pos'][i], b=tot_list['pos'][i+1], cpu=cpu, dims=n_input, verbose=verbose)

        negative_variants=[get_data(f_path+'neg.%d' %freq,a=tot_list[freq][i], b=tot_list[freq][i+1], cpu=cpu,verbose=verbose,dims=n_input) for freq in [0,5,10,15,20,25]]

        n_pos=np.vstack([tmp[0] for tmp in negative_variants])
        nx_train=np.vstack([tmp[1] for tmp in negative_variants])
        ny_train=np.vstack([tmp[2] for tmp in negative_variants])
        ntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
        ntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

        perm=np.random.permutation(len(nx_train))
        
        np.take(n_pos,perm,axis=0,out=n_pos)
        np.take(nx_train,perm,axis=0,out=nx_train)
        np.take(ny_train,perm,axis=0,out=ny_train)
        np.take(ntrain_allele,perm,axis=0,out=ntrain_allele)
        np.take(ntrain_ref,perm,axis=0,out=ntrain_ref)

        return (_,x_train,y_train,train_allele,train_ref), (n_pos,nx_train,ny_train,ntrain_allele,ntrain_ref)

    else:
        cpu=params['cpu']
        n_input=params['dims']
        
        if path:
            f_path=path
        else:
            f_path=params['test_path']
        
        _,vpx_train,vpy_train,vptrain_allele,vptrain_ref= get_data(f_path+'pos',cpu=cpu,dims=n_input)

        negative_variants=[get_data(f_path+'neg.%d' %freq,cpu=cpu,dims=n_input, verbose=verbose) for freq in [0,5,10,15,20,25]]

        vx_test=np.vstack([tmp[1] for tmp in negative_variants]+[vpx_train])[len(negative_variants[0][0])//2:]
        vy_test=np.vstack([tmp[2] for tmp in negative_variants]+[vpy_train])[len(negative_variants[0][0])//2:]
        vtest_allele=np.vstack([tmp[3] for tmp in negative_variants]+[vptrain_allele])[len(negative_variants[0][0])//2:]
        vtest_ref=np.vstack([tmp[4] for tmp in negative_variants]+[vptrain_ref])[len(negative_variants[0][0])//2:]
        
        return vx_test,vy_test,vtest_allele,vtest_ref            

def read_pileups_from_file(options):
    
    fname,n,mode,dims=options
    file= open(fname,'r')
    file.seek(n)
    mat=[]
    pos=[]
    ref=[]
    allele=[]
    gt=[]
    
    dp=[]
    freq=[]
    
    i=0
    if mode=='train':

        while i<1000:
            i+=1
            c=file.read(14)
            if not c:
                break
            pos.append(int(c[:11]))
            gt.append(int(c[11]))
            allele.append(int(c[12]))
            ref.append(int(c[13]))



            m=file.read(dims[0]*dims[1]*dims[2]*6)
            p_mat=np.array([int(m[6*i:6*i+6]) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))

            mat.append(p_mat)

        mat=np.array(mat)    
        pos=np.array(pos)
        ref=np.eye(4)[np.array(ref)].astype(np.int8)
        allele=np.eye(4)[np.array(allele)].astype(np.int8)
        gt=np.eye(2)[np.array(gt)].astype(np.int8)

        return (pos,mat,ref,allele,gt)

    else:
        while i<1000:
            i+=1
            c=file.read(20)
            if not c:
                break
            
            pos.append(int(c[:11]))
            dp.append(int(c[11:17]))
            freq.append(int(c[17:19]))
            ref.append(int(c[19]))



            m=file.read(dims[0]*dims[1]*dims[2]*6)
            p_mat=np.array([int(m[6*i:6*i+6]) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))

            mat.append(p_mat)

        mat=np.array(mat).astype(np.int16)    
        pos=np.array(pos)
        ref=np.eye(np.max(ref)+1)[np.array(ref)].astype(np.int8)
        ref=ref[:,:4]
        dp=np.array(dp)
        freq=np.array(freq)
        return (pos,mat,ref,dp,freq)



def pileup_image_from_mat(in_mat,size=(20,10)):
        t=np.copy(in_mat[:,:,-1])
        
        t=np.abs(t)
        t=t/np.max(t)
        const=np.min(t[t!=0])*0.5
        
        new=np.dstack([in_mat[:,:,:4],t,const*np.ones(t.shape)])
        p_mat=np.argmax(new,2)
        
        #deletion=purple,  1A=green,   2G=orange,    3T=red,    4C=blue
        color={5:np.array([[[0,0,0]]]),4:np.array([[[255, 255, 255]]]),0:np.array([[[0,   255, 0]]]),1:np.array([[[206, 132, 47]]]),2:np.array([[[255,   0, 0]]]),3:np.array([[[0,   0, 255]]])}
        data = np.zeros((p_mat.shape[0],p_mat.shape[1], 3)).astype(int)
        for j in range(6):
            new=(p_mat==j).astype(int)
            
            data+=np.repeat(new[:, :, np.newaxis], 3, axis=2)*color[j]
        
        plt.figure(figsize=(size[0],size[1]))
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
