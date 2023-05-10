from warnings import simplefilter 
simplefilter(action='ignore', category=FutureWarning)

import time, os, copy, argparse, subprocess, glob, re, datetime
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect_indel import *

if type(tf.contrib) != type(tf): tf.contrib._warning = None


config = tf.ConfigProto(device_count={"CPU": 32})
#config = tf.ConfigProto()
config.gpu_options.allow_growth = True


rev_gt_map={0:'hom-ref', 1:'hom-alt', 2:'het-ref', 3:'het-alt'}
rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

class Genome:
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.training_data = {}
        self.testing_data = None
        self.chunk_list=[]

    def add_data(self, pool, val):
        print('%s: Starting reading genome %s.' %(str(datetime.datetime.now()), self.name), flush=True)
        train_list=glob.glob(os.path.join(self.path,'train*'))
        chunk_list=[re.findall('pileups.(\d+)', x)[0] for x in train_list]
        
        print('\nProgress\n\nTotal: %s\nDone : ' %('.'*len(chunk_list)),end='', flush=True)
        
        for chunk, path in zip(chunk_list, train_list):
            self.chunk_list.append(chunk)
            x_train, train_gt= get_data(path,pool)

            self.training_data[chunk]=(x_train, train_gt)
                
            print('.',end='',flush=True)
        
        if val:
            self.testing_data = get_data(os.path.join(self.path,'test.pileups'),pool)
            
        print('\n\n%s: Finished reading genome %s.\n' %(str(datetime.datetime.now()), self.name), flush=True)
    
def train_indel_model(params):
    
    
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
            current_genome.add_data(pool, params['val'])
            
            genomes_list.append(current_genome)
            
    pool.close()  
    pool.join()

    print('\n%s: Finished reading pileups.' %str(datetime.datetime.now()),flush=True)
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    weights,biases,tensors=get_tensors([5,128,2],learning_rate)
    (x0, x1,x2,gt, accuracy, cost, optimizer, prob, rate)=tensors

    
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
                    
                    x_train, train_gt=genome.training_data[chunk]

                    for batch in range(len(x_train)//batch_size):
                        batch_x = x_train[batch*batch_size:min((batch+1)*batch_size, len(x_train))]
                        batch_gt = train_gt[batch*batch_size:min((batch+1)*batch_size, len(train_gt))]

                        opt, loss, batch_acc = sess.run([optimizer,cost,accuracy], feed_dict={x2: batch_x, gt:batch_gt, rate:0.5})

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
                        vx_test, vtest_gt= val_genome.testing_data
                        
                        v_loss,v_acc,total_test_data=0,0,0

                        for batch in range(int(np.ceil(len(vx_test)/(batch_size)))):
                            vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                            vbatch_gt = vtest_gt[batch*batch_size:min((batch+1)*batch_size,len(vtest_gt))] 

                            batch_loss,batch_acc= sess.run([cost, accuracy], feed_dict={x2: vbatch_x, gt:vbatch_gt, rate:0.0})
                            v_loss+=batch_loss*len(vbatch_x)
                            v_acc+=batch_acc
                            total_test_data+=len(vbatch_x)

                        print('Genome: %s    Validation Loss: %.4f    Validation Accuracy: %.4f' %(val_genome.name, v_loss/total_test_data, v_acc/total_test_data), flush=True)

                print(50*'*'+'\n', flush=True)
            saver.save(sess, save_path=os.path.join(params['model'],'model'),global_step=save_num)
            elapsed=time.time()-t
            save_num+=1
            
            print ('%s: Time Taken for Iteration %d: %.2f seconds\n' %(str(datetime.datetime.now()),k,elapsed), flush=True)
            
            

            
def get_data(fname,pool):
    t=time.time()
    l=os.stat(fname).st_size
    dims=(5,128,2)
    
    rec_size=12+dims[0]*dims[1]*dims[2]*4*3
    my_array=[(fname,x,'train',dims) for x in range(0,l,1000*rec_size)]
    
    results = pool.map(read_pileups_from_file, my_array)

    mat=np.vstack([res[0] for res in results])
    gt=np.vstack([res[1] for res in results])
    
    return mat,gt

def read_pileups_from_file(options):    
    fname,n,mode,dims=options
    file= open(fname,'r')
    file.seek(n)
    
    mat_0,mat_1,mat_2=[],[],[]
    pos=[]
    gt=[]
    
    i=0

    while i<1000:
        i+=1
        c=file.read(12)
        if not c:
            break             

        gt.append(int(c[11]))

        m=file.read(dims[0]*dims[1]*dims[2]*4)
        p_mat_0=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
        mat_0.append(p_mat_0)

        m=file.read(dims[0]*dims[1]*dims[2]*4)
        p_mat_1=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
        mat_1.append(p_mat_1)

        m=file.read(dims[0]*dims[1]*dims[2]*4)
        p_mat_2=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
        mat_2.append(p_mat_2)

    mat_0=norm(np.array(mat_0))
    mat_1=norm(np.array(mat_1))
    mat_2=norm(np.array(mat_2))

    mat=np.hstack([mat_0, mat_1, mat_2])

    gt=np.array(gt)
    gt=np.eye(4)[gt].astype(np.int8)

    return (mat,gt)


def norm(batch_x0):
    batch_x0=batch_x0.astype(np.float32)
    batch_x0[:,:,:,0]=batch_x0[:,:,:,0]/(np.sum(batch_x0[:,:,:,0],axis=1)[:,np.newaxis,:])-batch_x0[:,:,:,1]
    return batch_x0

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

    print('aasas',flush=True)
    os.makedirs(args.model, exist_ok=True)
    
    in_dict={'cpu':args.cpu,'rate':args.rate, 'iters':args.iterations, 'size':args.size,\
             'train_path':args.train, 'model':args.model, 'val':args.validation,'rt_path':args.retrain_path,\
            'sparse':args.sparse}
    
    train_indel_model(in_dict)

    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)