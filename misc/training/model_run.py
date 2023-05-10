from warnings import simplefilter 
simplefilter(action='ignore', category=FutureWarning)

import time, os, copy, argparse, subprocess, glob, re, datetime
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect import *
from utils import *
from scipy.sparse import csr_matrix

if type(tf.contrib) != type(tf): tf.contrib._warning = None


config = tf.ConfigProto(device_count={"CPU": 32})
#config = tf.ConfigProto()
config.gpu_options.allow_growth = True

num_to_base_map={0:'A',1:'G',2:'T',3:'C'}

class Genome:
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.training_data = {}
        self.testing_data = None
        self.chunk_list=[]

    def add_data(self, pool, sparse, val):
        print('%s: Starting reading genome %s.' %(str(datetime.datetime.now()), self.name), flush=True)
        train_list=glob.glob(os.path.join(self.path,'train*'))
        chunk_list=[re.findall('pileups.(\d+)', x)[0] for x in train_list]
        
        print('\nProgress\n\nTotal: %s\nDone : ' %('.'*len(chunk_list)),end='', flush=True)
        
        for chunk, path in zip(chunk_list, train_list):
            self.chunk_list.append(chunk)
            x_train,y_train,train_allele,train_ref= get_data(path,pool)

            if sparse:
                self.training_data[chunk]=(csr_matrix(x_train.reshape(-1)),y_train,train_allele,train_ref)
            else:
                self.training_data[chunk]=(x_train,y_train,train_allele,train_ref)
                
            print('.',end='',flush=True)
        
        if val:
            self.testing_data = get_data(os.path.join(self.path,'test.pileups'),pool)
            
        print('\n\n%s: Finished reading genome %s.\n' %(str(datetime.datetime.now()), self.name), flush=True)
    
def train_SNP_model(params):
    
    sparse=params['sparse']
    
    tf.reset_default_graph()
    cpu=params['cpu']
    params['val']=True
    dims=[5,41,5]
    n_input=dims
    
    pool = mp.Pool(processes=cpu)
    
    print('%s: Starting reading pileups.' %str(datetime.datetime.now()),flush=True)

    genomes_list=[]
    with open(params['train_path'],'r') as file:
        for line in file:
            x=line.rstrip('\n').split(',')
            
            current_genome=Genome(*x)
            current_genome.add_data(pool,sparse, params['val'])
            
            genomes_list.append(current_genome)
            
    pool.close()  
    pool.join()

    print('\n%s: Finished reading pileups.' %str(datetime.datetime.now()),flush=True)
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    weights,biases,tensors=get_tensors(n_input,learning_rate)
    (x,GT_label,A_label, G_label, T_label, C_label,GT_score, A_score, G_score, T_score, C_score, accuracy_GT, accuracy_A,  accuracy_G,  accuracy_T,  accuracy_C, prediction_accuracy_GT, prediction_accuracy_A,  prediction_accuracy_G,  prediction_accuracy_T,  prediction_accuracy_C, prediction_GT, prediction_A,  prediction_G,  prediction_T,  prediction_C, accuracy, cost, optimizer,  cost_GT, cost_A, cost_G, cost_T, cost_C,A_ref,G_ref,T_ref,C_ref,prob_GT,prob_A,prob_G,prob_T,prob_C,keep)=tensors

    init = tf.global_variables_initializer()
    saver = tf.train.Saver(max_to_keep=10000)
    
    
    with tf.Session(config=config)  as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        if params['rt_path']:
            print('%s: Retraining model:> %s' %(str(datetime.datetime.now()),params['rt_path']))
            saver.restore(sess, params['rt_path'])        

        stats,v_stats=[],[]
        
        count=0
        
        save_num=1        
        
        for k in range(training_iters):
            t=time.time()
            print('\n'+100*'-'+'\n',flush=True)
            print('%s: Starting epoch #: %d\n' %(str(datetime.datetime.now()),k),flush=True)
            
            for genome in genomes_list:
                print('%s: Training on genome %s \n' %(str(datetime.datetime.now()), genome.name), flush=True)
                
                chunk_list=genome.chunk_list
                print('Progress\n\nTotal: %s\nDone : ' %('.'*len(chunk_list)),end='', flush=True)
                
                for chunk in chunk_list:
                    training_loss, total_train_data, train_acc = 0, 0, 0
                    
                    x_train_sparse,y_train,train_allele,train_ref=genome.training_data[chunk]

                    if sparse:
                        x_train_flat=np.array(x_train_sparse.todense())
                        x_train=x_train_flat.reshape(x_train_flat.shape[1]//1025,5,41,5)
                    else:
                        x_train=x_train_sparse

                    for batch in range(len(x_train)//batch_size):
                        batch_x = x_train[batch*batch_size:min((batch+1)*batch_size, len(x_train))]

                        batch_y = y_train[batch*batch_size:min((batch+1)*batch_size, len(y_train))]

                        batch_ref = train_ref[batch*batch_size :min((batch+1)*batch_size, len(train_ref))]

                        batch_allele = train_allele[batch*batch_size :min((batch+1)*batch_size, len(train_allele))]

                        opt,loss,batch_acc = sess.run([optimizer,cost,accuracy], feed_dict={x: batch_x, GT_label:batch_y, A_label:np.eye(2)[batch_allele[:,0]], G_label:np.eye(2)[batch_allele[:,1]], T_label:np.eye(2)[batch_allele[:,2]], C_label:np.eye(2)[batch_allele[:,3]] , A_ref:batch_ref[:,0][:,np.newaxis], G_ref:batch_ref[:,1][:,np.newaxis], T_ref:batch_ref[:,2][:,np.newaxis], C_ref:batch_ref[:,3][:,np.newaxis], keep:0.5})

                        training_loss+=loss*len(batch_x)
                        total_train_data+=len(batch_x)
                        train_acc+=batch_acc

                    

                    print('.',end='',flush=True)

                    training_loss=training_loss/total_train_data
                    train_acc=train_acc/total_train_data
                    
                print('\n\nGenome: %s      Training Loss: %.4f      Training Accuracy: %.4f\n' %(genome.name, training_loss, train_acc), flush=True)
                
                if params['val']:
                    print(50*'*'+'\n', flush=True)
                    print('%s: Performing Validation\n' %str(datetime.datetime.now()), flush=True)
                    
                    for val_genome in genomes_list:
                        vx_test, vy_test, vtest_allele, vtest_ref = val_genome.testing_data
                        
                        v_loss,v_acc,A_acc,G_acc,T_acc,C_acc,GT_acc,v_loss,v_acc=0,0,0,0,0,0,0,0,0

                        for batch in range(int(np.ceil(len(vx_test)/(batch_size)))):
                            vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                            vbatch_y = vy_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))] 
                            vbatch_allele = vtest_allele[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                            vbatch_ref = vtest_ref[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]

                            batch_loss,batch_acc, batch_GT_acc, batch_A_acc, batch_G_acc, batch_T_acc, batch_C_acc,batch_prediction_GT, batch_prediction_A,  batch_prediction_G,  batch_prediction_T,  batch_prediction_C,bGT_score,bA_score,bG_score,bT_score,bC_score = sess.run([cost, accuracy, accuracy_GT, accuracy_A,  accuracy_G,  accuracy_T,  accuracy_C,prediction_GT, prediction_A,  prediction_G,  prediction_T,  prediction_C,GT_score,A_score,G_score,T_score,C_score], feed_dict={x: vbatch_x,GT_label:vbatch_y, A_label:np.eye(2)[vbatch_allele[:,0]], G_label:np.eye(2)[vbatch_allele[:,1]], T_label:np.eye(2)[vbatch_allele[:,2]], C_label:np.eye(2)[vbatch_allele[:,3]],  A_ref:vbatch_ref[:,0][:,np.newaxis], G_ref:vbatch_ref[:,1][:,np.newaxis], T_ref:vbatch_ref[:,2][:,np.newaxis], C_ref:vbatch_ref[:,3][:,np.newaxis], keep:1.0})

                            v_loss+=batch_loss*len(vbatch_x)
                            v_acc+=batch_acc
                            A_acc+=batch_A_acc
                            G_acc+=batch_G_acc
                            T_acc+=batch_T_acc
                            C_acc+=batch_C_acc
                            GT_acc+=batch_GT_acc

                        print('Genome: %s    Validation Loss: %.4f    Validation Accuracy: %.4f' %(val_genome.name, v_loss/len(vx_test), v_acc/len(vx_test)), flush=True)

                        print('Genome: %s    GT_acc=%.4f    A_acc=%.4f    G_acc=%.4f    T_acc=%.4f    C_acc=%.4f\n' %(val_genome.name, GT_acc/len(vx_test), A_acc/len(vx_test), G_acc/len(vx_test), T_acc/len(vx_test), C_acc/len(vx_test)), flush=True)

                print(50*'*'+'\n', flush=True)
            saver.save(sess, save_path=os.path.join(params['model'],'model'),global_step=save_num)
            elapsed=time.time()-t
            save_num+=1
            
            print ('%s: Time Taken for Iteration %d: %.2f seconds\n' %(str(datetime.datetime.now()),k,elapsed), flush=True)
            
            

            
def get_data(fname,pool,a=None, b=None,dims=(5,41,5)):
    t=time.time()
    l=os.stat(fname).st_size
    
    rec_size=15+dims[0]*dims[1]*dims[2]*6
    if a!=None and b!=None:
        my_array=[(fname,x,dims) for x in range(a,b,1000*rec_size)]
    else:
        my_array=[(fname,x,dims) for x in range(0,l,1000*rec_size)]

    results = pool.map(read_pileups_from_file, my_array)
    
    pos=np.vstack([res[0][:,np.newaxis] for res in results])
    mat=np.vstack([res[1] for res in results])
    

    ref=np.vstack([res[2] for res in results])
    allele, gt=None, None

    allele=np.vstack([res[3] for res in results])
    
    gt=np.vstack([res[4] for res in results])
   

    return mat,gt,allele,ref
    
def read_pileups_from_file(options):
    
    fname,n,dims=options
    file= open(fname,'r')
    file.seek(n)
    mat=[]
    pos=[]
    ref=[]
    allele1,allele2=[],[]
    gt=[]
    
    dp=[]
    freq=[]
    
    i=0
    while i<1000:
        i+=1
        c=file.read(15)
        if not c:
            break
        pos.append(int(c[:11]))
        gt.append(int(c[11]))
        allele1.append(int(c[12]))
        allele2.append(int(c[13]))
        ref.append(int(c[14]))



        m=file.read(dims[0]*dims[1]*dims[2]*6)
        p_mat=np.array([int(m[6*i:6*i+6]) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))

        mat.append(p_mat)

    mat=np.array(mat)    
    pos=np.array(pos)
    ref=np.eye(4)[np.array(ref)].astype(np.int8)
    allele1=np.eye(4)[np.array(allele1)].astype(np.int8)
    allele2=np.eye(4)[np.array(allele2)].astype(np.int8)
    
    allele=allele1+allele2
    
    allele[allele==2]=1
    
    gt=np.eye(2)[np.array(gt)].astype(np.int8)

    return (pos,mat.astype(np.int16),ref,allele,gt)

if __name__ == '__main__':
    t=time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rate", help="Learning rate",type=float)
    parser.add_argument("-i", "--iterations", help="Training iterations",type=int)
    parser.add_argument("-s", "--size", help="Batch size",type=int)
    parser.add_argument("-train", "--train", help="Train path")
    parser.add_argument("-model", "--model", help="Model output path")
    parser.add_argument("-cpu", "--cpu", help="CPUs",type=int)
    parser.add_argument("-rt_path", "--retrain_path", help="Retrain saved model",type=str)
    parser.add_argument('-sparse','--sparse', help='Stores features as sparse matrices', default=False, action='store_true')
    parser.add_argument('-val','--validation', help='Perform validation', default=True, action='store_false')
    
    args = parser.parse_args()

    
    os.makedirs(args.model, exist_ok=True)
    
    in_dict={'cpu':args.cpu,'rate':args.rate, 'iters':args.iterations, 'size':args.size,\
             'train_path':args.train, 'model':args.model, 'val':args.validation,'rt_path':args.retrain_path,\
            'sparse':args.sparse}
    
    train_SNP_model(in_dict)

    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)