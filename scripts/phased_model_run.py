import time,os,copy,argparse,subprocess,psutil
import pandas as pd
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from phased_model_architect import *

config = tf.ConfigProto(device_count={"CPU": 64})
#config = tf.ConfigProto()
config.gpu_options.allow_growth = True


def genotype_caller_skinny(params,input_type='path',data=None,attempt=0,neg_part='neg.combined'):
    tf.reset_default_graph()
    
    false_mltplr=3
    
    cpu=params['cpu']
    n_input=params['dims']
    dims=n_input
    chrom_list=list(range(2,23)) #params['chrom'].split(':') 
    #chrom_list=list(range(int(chrom_list[0]),int(chrom_list[1])+1))
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    weights,biases,tensors=get_tensors(n_input,learning_rate)
    (x,allele,ref,fc_layer,cost,optimizer,keep,correct_prediction,prediction, accuracy)=tensors
    
    if params['val']:
        val_list=[]
        for v_path in params['test_path'].split(':'):
            vx_test, vtest_allele, vtest_ref=get_data_20plus(params,v_path)
            val_list.append((v_path,(vx_test, vtest_allele, vtest_ref)))
        
    
    init = tf.global_variables_initializer()
    saver = tf.train.Saver(max_to_keep=100)
    
    rec_size=13+dims[0]*dims[1]*dims[2]*3
    
    
    n_size=1
    with tf.Session(config=config)  as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        if params['retrain']:
            saver.restore(sess, params['rt_path'])        

        stats,v_stats=[],[]
        
        count=0
        
        save_num=1
        t=time.time()
        
        iter_ratio=params['ratio'] if params['ratio'] else 20
        
        iter_steps=max(training_iters//iter_ratio,1)
        iters=min(iter_ratio,training_iters)
        chrom_data={}
        
        print('Starting reading pileups',flush=True)
        for chrom in chrom_list:
                f_path=os.path.join(params['train_path'],'chr%d/chr%d.pileups.' %(chrom,chrom))
                                    
                _,x_train,train_allele,train_ref= get_data_ps(f_path+'pos',cpu=cpu, dims=n_input)
                _,nx_train,ntrain_allele,ntrain_ref=get_data_ps(f_path+neg_part, cpu=cpu, dims=n_input)
                chrom_data['chr%d'%chrom]=((x_train,train_allele,train_ref),(nx_train,ntrain_allele,ntrain_ref))
                print('chromosome %d done | '%chrom,end='',flush=True)

        print('Finished reading pileups',flush=True)
        
        for k in range(iter_steps):
            
            for chrom in chrom_list:
                print('Training on chrom %d ' %(chrom),end='',flush=True)
                (x_train,train_allele,train_ref), (nx_train,ntrain_allele,ntrain_ref)=chrom_data['chr%d'%chrom]


                n_start=-false_mltplr*len(x_train)
                n_end=0
                for i in range(iters):


                    n_start+=false_mltplr*len(x_train)
                    n_end=n_start+false_mltplr*len(x_train)

                    if n_end>len(nx_train):
                        batch_nx_train= np.vstack([nx_train[n_start:,:,:,:],nx_train[:n_end-len(nx_train),:,:,:]]) 

                        batch_ntrain_allele=np.vstack([ntrain_allele[n_start:,:],ntrain_allele[:n_end-len(nx_train),:]])
                        batch_ntrain_ref = np.vstack([ntrain_ref[n_start:,:],ntrain_ref[:n_end-len(nx_train),:]])


                        n_start=n_end-len(nx_train)
                        n_end=n_start+false_mltplr*len(x_train)
                    else:    
                        batch_nx_train,batch_ntrain_allele, batch_ntrain_ref = \
                        nx_train[n_start:n_end,:,:,:],\
                        ntrain_allele[n_start:n_end,:], ntrain_ref[n_start:n_end,:]


                    training_loss=0
                    total_train_data=0

                    total_false_size=int(false_mltplr*batch_size)

                    for batch in range(int(np.ceil(len(x_train)/batch_size))):
                        batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))],\
                                  batch_nx_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx_train))]])


                        batch_allele = np.vstack([train_allele[batch*batch_size :min((batch+1)*batch_size,\
                                                   len(train_allele))], batch_ntrain_allele[total_false_size*batch : \
                                                   min(total_false_size*(batch+1), len(batch_ntrain_allele))]])

                        batch_ref = np.vstack([train_ref[batch*batch_size :min((batch+1)*batch_size,\
                                                len(train_ref))], batch_ntrain_ref[ total_false_size*batch:\
                                                min(total_false_size*(batch+1), len(batch_ntrain_ref))]])

                        opt,loss = sess.run([optimizer,cost], feed_dict={x: batch_x, ref:batch_ref, allele:batch_allele, keep:0.5})
                        training_loss+=loss*len(batch_x)
                        total_train_data+=len(batch_x)

                    training_loss=training_loss/total_train_data

                    print('.',end='',flush=True)

                if params['val'] and (k<3 or chrom==22):

                    for val in val_list:
                        print('\n')
                        print(30*'-')
                        print(val[0])
                        vx_test, vtest_allele, vtest_ref=val[1]


                        v_loss,v_acc=0,0

                        v_score,v_pred=[],[]

                        for batch in range(int(np.ceil(len(vx_test)/(batch_size)))):
                            vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                            vbatch_allele = vtest_allele[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                            vbatch_ref = vtest_ref[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]

                            batch_loss,batch_acc=sess.run([cost,accuracy], feed_dict={x: vbatch_x, ref:vbatch_ref, allele:vbatch_allele, keep:1.0})

                            v_loss+=batch_loss*len(vbatch_x)
                            v_acc+=batch_acc
                     

                        print('train loss= %.4f' %loss)
                        print('valid loss= %.4f\n' %(v_loss/len(vx_test)))
                        print('valid accuracy= %.4f' %(v_acc/len(vx_test)))
                        
                        print(100*'.'+'\n')


            saver.save(sess, save_path=params['model'],global_step=save_num)
            elapsed=time.time()-t
            
            print ('Time Taken for Iteration %d-%d: %.2f seconds\n'\
                   %((save_num-1)*iters,save_num*iters,elapsed), flush=True)
            
            save_num+=1
            t=time.time()
            
        
        

def test_model(params,suffix='',prob_save=False):
    model_path,test_path,n_input,chrom,vcf_path= params['model'], params['test_path'],params['dims'],params['chrom'],params['vcf_path']
    
    cpu=params['cpu']
    tf.reset_default_graph()
    
    dims=n_input[:]
    params['window']=None

    weights,biases,tensors=get_tensors(n_input,0.0)
    (x,allele,ref,fc_layer,cost,optimizer,keep,probability, accuracy)=tensors

    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    
    rev_mapping={0:'A',1:'G',2:'T',3:'C'}
    ts_tv={0:5,1:7,2:2,3:3}
    gt_map={0:1,1:0}
    rec_size=12+dims[0]*dims[1]*dims[2]*6
    batch_size=1000
    total=[]
    
    neg_file=open('%s.%s.stats' %(vcf_path,params['phase']),'w')
    vcf_path=vcf_path+suffix
    
    with open(vcf_path,'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        

        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GP,Number=1,Type=Integer,Description="Genotype Probability">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
        ttt=[]
        
        chnk=1
        tot_list={}
        f_path=params['test_path']
        
        statinfo = os.stat(f_path)
        sz=statinfo.st_size
        tmp_sz=list(range(0,sz,rec_size*(sz//(chnk*rec_size))))
        tmp_sz=tmp_sz[:chnk]
        tmp_sz=tmp_sz+[sz] if tmp_sz[-1]!=sz else tmp_sz
        
        threshold=0.5
        for i in range(len(tmp_sz)-1):
            pos,x_test,test_ref,test_dp= get_data_ps(f_path,a=tmp_sz[i], b=tmp_sz[i+1],dims=n_input,cpu=cpu,mode='test')
            
            z=np.sum(np.sum(np.abs(x_test[1:,:,:4]),axis=2),axis=0)
            z=z[:,np.newaxis]
            z[z==0]=1
            x_test[1:,:,:4]=x_test[1:,:,:4]/z
            
            
            for batch in range(int(np.ceil(len(x_test)/batch_size))):
                        batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                        batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]

                        batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]
                        batch_dp = test_dp[batch*batch_size:min((batch+1)*batch_size, len(test_dp))]

                        batch_probs= sess.run(probability, feed_dict={x: batch_x,keep:1.0})
                                                
                        batch_pred=np.argmax(batch_probs,1)
                        batch_max_prob=np.max(batch_probs,1)
                        batch_ref=np.argmax(batch_ref,1)

                        for j in range(len(batch_ref)):
                            
                            if batch_pred[j]!=batch_ref[j]:
                                s='%s\t%d\t.\t%s\t%s\t%d\t%s\t.\tGT:GP:DP\t%s:%d:%d\n' %(chrom, batch_pos[j], rev_mapping[batch_ref[j]], rev_mapping[batch_pred[j]], 0,'PASS','1', int(min(99,-10*np.log10(1e-10+ 1-batch_max_prob[j]))), batch_dp[j])
                               
                                f.write(s)

                        neg_file.write('%d,%.4f,%.4f,%.4f,%.4f,%d\n' %(batch_pos[j], batch_probs[j,0], batch_probs[j,1], batch_probs[j,2], batch_probs[j,3],batch_dp[j]))
                                                                                                   
                                                                                                                   
                       
    neg_file.close()
    return vcf_path
'''
def test_with_hap(params):
    count=1
    change=1
    init_vcf_path=test_model(params,suffix='.initial')
    tmp_vcf_path=init_vcf_path
    
    while True:
        tmp_nbr_params={}
        nbr_path=get_neighbors.generate(tmp_nbr_params,mode='redo',count)
        
        tmp_cnd_params={}
        new_cand_path=gcp.generate()
        
        tmp_test_params={}
        tmp_vcf_path
        
        change=0
        
        if change<0.1:
            #name change
            break
    
'''        

def make_allele(ylabel,allele,reference):
    new_allele=allele.copy()
    new_allele[ylabel[:,1]==1]=new_allele[ylabel[:,1]==1]+reference[ylabel[:,1]==1]
    return new_allele

def get_data_ps(fname,a=None, b=None,dims=(32,33,5), cpu=4,mode='train',verbose=False):
    t=time.time()
    l=os.stat(fname).st_size
    
    if mode=='train':
        rec_size=13+dims[0]*dims[1]*dims[2]*3
        if a!=None and b!=None:
            my_array=[(fname,x,mode,dims) for x in range(a,b,1000*rec_size)]
        else:
            my_array=[(fname,x,mode,dims) for x in range(0,l,1000*rec_size)]
    else:
        rec_size=16+dims[0]*dims[1]*dims[2]*3
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
    allele=None
    if mode=='train':
        allele=np.vstack([res[3] for res in results])
        return pos,mat,allele,ref

    else:
        dp=np.vstack([res[3] for res in results])
        return pos,mat,ref,dp
    
def get_data_20plus(params,f_path):
    cpu=params['cpu']
    n_input=params['dims']
    
    _,vpx_train,vptrain_allele,vptrain_ref= get_data_ps(f_path+'pos',cpu=cpu,dims=n_input)
    _,vnx_test,vntest_allele,vntest_ref=get_data_ps(f_path+'neg.combined',cpu=cpu,dims=n_input)
    vx_test,vtest_allele,vtest_ref =np.vstack([vpx_train,vnx_test]), np.vstack([vptrain_allele,vntest_allele]), np.vstack([vptrain_ref,vntest_ref])
    return (vx_test,vtest_allele,vtest_ref)

    
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
    i=0
    if mode=='train':

        while i<1000:
            i+=1
            c=file.read(13)
            if not c:
                break
            try:
                pos.append(int(c[:11]))
                allele.append(int(c[11]))
                ref.append(int(c[12]))

            except ValueError:
                print(n,i,fname)
                print(dims)
                pos.append(int(c[:11]))
                allele.append(int(c[11]))
                ref.append(int(c[12]))
                


            m=file.read(dims[0]*dims[1]*dims[2]*3)
            p_mat=np.array([int(m[3*i:3*i+3]) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0], dims[1], dims[2]))

            mat.append(p_mat)

        mat=np.array(mat)    
        pos=np.array(pos)
        ref=np.eye(4)[np.array(ref)].astype(np.int8)
        allele=np.eye(4)[np.array(allele)].astype(np.int8)

        return (pos,mat.astype(np.int8),ref,allele)

    else:
        while i<1000:
            i+=1
            c=file.read(16)
            if not c:
                break
            
            try:
                pos.append(int(c[:11]))
                ref.append(int(c[11]))
                dp.append(int(c[12:]))



                m=file.read(dims[0]*dims[1]*dims[2]*3)
                p_mat=np.array([int(m[3*i:3*i+3]) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))

                mat.append(p_mat)
            except ValueError:
                print(i,n)
        mat=np.array(mat)    
        pos=np.array(pos)
        ref=np.eye(max(4,np.max(ref)+1))[np.array(ref)].astype(np.int8)
        ref=ref[:,:4]
        dp=np.array(dp)[:,np.newaxis]
        
        return (pos,mat.astype(np.int8),ref,dp.astype(np.int16))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rate", help="Learning rate",type=float)
    parser.add_argument("-i", "--iterations", help="Training iterations",type=int)
    parser.add_argument("-s", "--size", help="Batch size",type=int)
    parser.add_argument("-train", "--train", help="Train path")
    parser.add_argument("-test", "--test", help="Test path")
    parser.add_argument("-model", "--model", help="Model output path")
    parser.add_argument("-m", "--mode", help="Mode")
    parser.add_argument("-dim", "--dimensions", help="Input dimensions")
    parser.add_argument("-vcf", "--vcf", help="VCF output path")
    parser.add_argument("-chrom", "--chrom", help="Chromosome")
    parser.add_argument("-cpu", "--cpu", help="CPUs",type=int)
    parser.add_argument("-val", "--validation", help="Validation",type=int)
    parser.add_argument("-rt", "--retrain", help="Retrain saved model",type=int)
    parser.add_argument("-w", "--window", help="Window size around site",type=int)
    parser.add_argument("-ratio", "--ratio", help="iterations per batch",type=int)
    parser.add_argument("-neg", "--neg_part", help="Negative Part")
    parser.add_argument("-wdir", "--workdir", help="Working Directory")
    parser.add_argument("-rt_path", "--rt_path", help="re-train directory")
    parser.add_argument("-phase", "--phase", help="Phase set")
    
    
    args = parser.parse_args()
    input_dims=[int(x) for x in args.dimensions.split(':')]
    t=time.time()
    
    if args.mode=='train':
        in_dict={'cpu':args.cpu,'rate':args.rate, 'iters':args.iterations, 'size':args.size,'dims':input_dims,'chrom':args.chrom,\
                 'train_path':args.train, 'test_path':args.test, 'model':args.model, 'val':args.validation,'retrain':args.retrain,\
                'window':args.window,'ratio':args.ratio,'rt_path':args.rt_path}
        genotype_caller_skinny(in_dict,neg_part=args.neg_part)
    
    else:
        in_dict={'cpu':args.cpu,'dims':input_dims,'test_path':args.test,'model':args.model, 'chrom':args.chrom, 'vcf_path':args.vcf, 'workdir':args.workdir,'phase':args.phase}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
