import time,os,copy,argparse,utils,gzip
import pandas as pd
import numpy as np
import multiprocessing as mp
import generate_pileups as gcp
import tensorflow as tf

def get_data(fname,cpu=4,mode='train'):
    t=time.time()
    file= open(fname,'r')
    l=len(file.read())
    if mode=='train':
        my_array=[(fname,x,mode) for x in range(0,l,1000*7406)]
    else:
        my_array=[(fname,x,mode) for x in range(0,l,1000*7404)]

    pool = mp.Pool(processes=cpu)
    results = pool.map(utils.read_pileups_from_file, my_array)
    elapsed=time.time()-t
    print ('I/O Time Elapsed: %.2f seconds' %elapsed,flush=True)
    
    pos=np.vstack([res[0][:,np.newaxis] for res in results])
    mat=np.vstack([res[1] for res in results])
    ref=np.vstack([res[2] for res in results])
    
    allele,gt=None,None
    if mode=='train':
        allele=np.vstack([res[3] for res in results])
        gt=np.vstack([res[4] for res in results])
    
    return pos,mat,gt,allele,ref
        
        
def conv2d(x, W, b, strides=1):
    x = tf.nn.conv2d(x, W, strides=[1, strides, strides, 1], padding='VALID')
    x = tf.nn.bias_add(x, b)
    return tf.nn.selu(x)




def conv_net(x,ref, weights, biases,keep):  
    conv1 = conv2d(x, weights['wc1'], biases['bc1'],2)

    conv2 = conv2d(conv1, weights['wc2'], biases['bc2'],2)

    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],2)

    flat_nn = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    
    drop_out_flat =flat_nn# tf.nn.dropout(flat_nn, keep)
    
    fc1 = tf.add(tf.matmul(drop_out_flat, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.selu(fc1)
    
    drop_out_fc1=tf.nn.dropout(fc1, keep)
    
    fc2 = tf.add(tf.matmul(drop_out_fc1, weights['wd-gtp']), biases['bd-gtp'])
    fc2 = tf.nn.selu(fc2)
    out_gt= tf.add(tf.matmul(fc2, weights['out-gtp']), biases['out-gtp'])
    
    fa = tf.add(tf.matmul(tf.concat([fc2,drop_out_fc1,ref],1), weights['w-tr']), biases['w-tr'])
    fa = tf.nn.selu(fa)
    out_allele = tf.add(tf.matmul(fa, weights['out-all']), biases['out-all'])

    return out_gt,out_allele

def get_tensors(n_input,learning_rate=0):
    h_in,w_in,depth=n_input
    for i in range(3):
        h_in=int(np.ceil(float(h_in-3+1)/float(2)))
        w_in=int(np.ceil(float(w_in-3+1)/float(2)))
        
    weights = {
        'wc1': tf.get_variable('W0', shape=(3,3,depth,32), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc2': tf.get_variable('W1', shape=(3,3,32,64), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(3,3,64,128), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(h_in*w_in*128,64), initializer=tf.contrib.layers.xavier_initializer()), 
        'w-tr': tf.get_variable('W4', shape=(100,32), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('W5', shape=(32,4), initializer=tf.contrib.layers.xavier_initializer()),
        'wd-gtp': tf.get_variable('W6', shape=(64,32), initializer=tf.contrib.layers.xavier_initializer()), 
        'out-gtp': tf.get_variable('W7', shape=(32,2), initializer=tf.contrib.layers.xavier_initializer()), 

    
    }
    
    biases = {
        'bc1': tf.get_variable('B0', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(64), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(64), initializer=tf.contrib.layers.xavier_initializer()),
        'w-tr': tf.get_variable('B4', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('B5', shape=(4), initializer=tf.contrib.layers.xavier_initializer()),
        'bd-gtp': tf.get_variable('B6', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'out-gtp': tf.get_variable('B7', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),

    
    }

    x,y = tf.placeholder("float", [None]+n_input), tf.placeholder("float", [None, 2])
    allele,ref=tf.placeholder("float", [None, 4]),tf.placeholder("float", [None, 4])
    keep = tf.placeholder(tf.float32)
    
    pred,fc_layer = conv_net(x,ref,weights, biases,keep)
    
    gt_likelihood=tf.nn.softmax(logits=pred)
    allele_likelihood=tf.nn.softmax(logits=fc_layer)
    
    cost_gt=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=pred,\
    labels=y))
    cost_allele=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=fc_layer, labels=allele))
    gamma=0.000
    trade_off=1.0
    #reg_loss=tf.add_n([tf.nn.l2_loss(t) for t in weights.values()])
                      #+[tf.constant(0.005)*tf.nn.l2_loss(weights['w-tr'])+tf.constant(0.005)*tf.nn.l2_loss(weights['out-all'])])
                      
    cost = cost_allele+cost_gt#+tf.constant(gamma)*tf.reduce_mean(reg_loss)
    
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    
    correct_prediction_gt = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
    accuracy_gt = tf.reduce_mean(tf.cast(correct_prediction_gt, tf.float32))
    
    correct_prediction_allele = tf.equal(tf.argmax(fc_layer, 1), tf.argmax(allele, 1))
    accuracy_allele = tf.reduce_mean(tf.cast(correct_prediction_allele, tf.float32))
    
    t1=(x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)
    t2=(correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)
    return weights,biases,t1,t2
    


    
def genotype_caller(params,input_type='path',data=None):
    tf.reset_default_graph()
    cpu=params['cpu']
    n_input=params['dims']
    if input_type=='path':
        _,x_train,y_train,train_allele,train_ref= get_data(params['train_path']+'pos',cpu=cpu)
        _,nx_train,ny_train,ntrain_allele,ntrain_ref= get_data(params['train_path']+'neg',cpu=cpu)
        
        
        _,vpx_test,vpy_test,vptest_allele,vptest_ref=get_data(params['test_path']+'pos',cpu=cpu)
        _,vnx_test,vny_test,vntest_allele,vntest_ref=get_data(params['test_path']+'neg',cpu=cpu)

        test_size=len(vpx_test)//10
        vx_test,vy_test,vtest_allele,vtest_ref=np.vstack([vpx_test,vnx_test[:test_size]]),\
        np.vstack([vpy_test,vny_test[:test_size]]),np.vstack([vptest_allele,vntest_allele[:test_size]]),\
        np.vstack([vptest_ref,vntest_ref[:test_size]])      
        print('train data received')
    else:
        _,x_train,y_train,train_allele,train_ref=data['train_data'][0]
        _,nx_train,ny_train,ntrain_allele,ntrain_ref= data['train_data'][1]
        
        _,vpx_test,vpy_test,vptest_allele,vptest_ref=data['test_data'][0]
        _,vnx_test,vny_test,vntest_allele,vntest_ref=data['test_data'][1]
        
        test_size=len(vpx_test)//10
        vx_test,vy_test,vtest_allele,vtest_ref=np.vstack([vpx_test[:test_size],vnx_test[:test_size]]),\
        np.vstack([vpy_test[:test_size],vny_test[:test_size]]),np.vstack([vptest_allele[:test_size],vntest_allele[:test_size]]),\
        np.vstack([vptest_ref[:test_size],vntest_ref[:test_size]])      
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    n_input=[i for i in x_train.shape[1:]]

    weights,biases,t1,t2=get_tensors(n_input,learning_rate)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2

    init = tf.global_variables_initializer()
    saver = tf.train.Saver(max_to_keep=1)

    
    n_size=1
    with tf.Session() as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        stats=[]
        v_stats=[]
        print('starting training',flush=True)
        count=0
        n_start=-len(x_train)
        n_end=0
        best_error=np.inf
        loss_list=[np.inf]
        p=5
        j=0
        c=0
        for i in range(training_iters):
            n_start+=len(x_train)
            n_end=n_start+len(nx_train)
            if n_end>len(x_train):
                n_start=0
                n_end=n_start+len(nx_train)
            batch_nx_train,batch_ny_train,batch_ntrain_allele,batch_ntrain_ref=nx_train[n_start:n_end,:,:,:],\
                ny_train[n_start:n_end,:],ntrain_allele[n_start:n_end,:],ntrain_ref[n_start:n_end,:]
            for batch in range(len(x_train)//batch_size):
                
                batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))],batch_nx_train[batch*batch_size:min((batch+1)*batch_size,len(batch_nx_train))]])
                batch_y = np.vstack([y_train[batch*batch_size:min((batch+1)*batch_size,len(y_train))],batch_ny_train[batch*batch_size:min((batch+1)*batch_size,len(batch_ny_train))]])    
                batch_ref = np.vstack([train_ref[batch*batch_size:min((batch+1)*batch_size,len(train_ref))],batch_ntrain_ref[batch*batch_size:min((batch+1)*batch_size,len(batch_ntrain_ref))]])
                batch_allele = np.vstack([train_allele[batch*batch_size:min((batch+1)*batch_size,len(train_allele))],batch_ntrain_allele[batch*batch_size:min((batch+1)*batch_size,len(batch_ntrain_allele))]])
                # Run optimization op (backprop).
                    # Calculate batch loss and accuracy
                opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele,keep:0.5})
                loss, acc_gt,acc_allele,fc_layer_batch,score_batch = sess.run([cost, accuracy_gt,accuracy_allele,fc_layer,pred], feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele,keep:1.0})
                
                stats.append([loss, acc_gt,acc_allele])
                
            all_pred,gt_pred=np.argmax(fc_layer_batch,axis=1),np.argmax(score_batch,axis=1)
            train_mat=np.hstack([all_pred[:,np.newaxis],np.argmax(batch_allele,axis=1)[:,np.newaxis],\
                                           np.argmax(batch_ref,axis=1)[:,np.newaxis], gt_pred[:,np.newaxis],\
                                           np.argmax(batch_y,axis=1)[:,np.newaxis]])
            
            train_acc=sum((train_mat[:,0]==train_mat[:,1])&(train_mat[:,4]==train_mat[:,3]))/len(train_mat)
            
            loss_list.append(loss)
            loss_diff=max(loss_list[-p:])-min(loss_list[-p:])
            print('train loss= %.4f   train loss diff= %.4f\n' %(loss,loss_diff))
            
            
            
            test_mat=[]
            for batch in range(len(vx_test)//batch_size):
                vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                vbatch_y = vy_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))] 
                vbatch_ref = vtest_ref[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                vbatch_allele = vtest_allele[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                
                
                fc_layer_batch,score_batch,v_loss = sess.run([fc_layer,pred,cost], feed_dict={x: vbatch_x,y: vbatch_y,ref:vbatch_ref,allele:vbatch_allele,keep:1.0})
                
                all_pred,gt_pred=np.argmax(fc_layer_batch,axis=1),np.argmax(score_batch,axis=1)

                test_mat.append(np.hstack([all_pred[:,np.newaxis],np.argmax(vbatch_allele,axis=1)[:,np.newaxis],\
                                           np.argmax(vbatch_ref,axis=1)[:,np.newaxis], gt_pred[:,np.newaxis],\
                                           np.argmax(vbatch_y,axis=1)[:,np.newaxis]]))
            
            test_mat=np.vstack(test_mat)
            #test_mat=test_mat[(test_mat[:,0]!=test_mat[:,1]) | (test_mat[:,2]!=test_mat[:,1])]
            v_acc=sum((test_mat[:,0]==test_mat[:,1])&(test_mat[:,4]==test_mat[:,3]))/len(test_mat)
            v_gt_acc=sum(test_mat[:,4]==test_mat[:,3])/len(test_mat)
            v_allele_acc=sum(test_mat[:,0]==test_mat[:,1])/len(test_mat)
            v_error=1-v_acc
            v_stats.append(v_error)
            
            
            print('train accuracy= %.4f           valid accuracy= %.4f' %(train_acc,v_acc))
            print('train GT accuracy= %.4f        valid GT accuracy= %.4f' %(acc_gt,v_gt_acc))
            print('train Allele accuracy= %.4f    valid Allele accuracy= %.4f' %(acc_allele,v_allele_acc))
            print(50*'.')
            print('\n')
            if v_error<best_error:
                best_error=v_error
                j=0
                saver.save(sess, save_path=params['model'])
                c=i
            elif loss_diff<1e-4:
                j+=1
            if j>=p:
                break
    print('stopping at iteration number: %d' %c,flush=True)            
    stats=np.array(stats)
    loss_list=np.array(loss_list)
    np.save('stats_full',stats)
    v_stats=np.array(v_stats)
    np.savez('v_stats_full',v_stats=v_stats,loss_list=loss_list)
    return stats,v_stats,loss_list

def genotype_caller_skinny(params,input_type='path',data=None):
    tf.reset_default_graph()
    cpu=params['cpu']
    n_input=params['dims']
    '''if input_type=='path':
        _,x_train,y_train,train_allele,train_ref= get_data(params['train_path']+'pos',cpu=cpu)
        _,nx_train,ny_train,ntrain_allele,ntrain_ref= get_data(params['train_path']+'neg',cpu=cpu)

        print('train data received',flush=True)
    else:
        _,x_train,y_train,train_allele,train_ref=data['train_data'][0]
        _,nx_train,ny_train,ntrain_allele,ntrain_ref= data['train_data'][1]'''
      
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']


    weights,biases,t1,t2=get_tensors(n_input,learning_rate)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2

    init = tf.global_variables_initializer()
    saver = tf.train.Saver(max_to_keep=10)

    
    n_size=1
    with tf.Session() as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        stats=[]
        v_stats=[]
        print('starting training',flush=True)
        count=0
        
        save_num=1
        t=time.time()
        for k in range(training_iters//10):
            
            for chrom in range(2,21):
                f_path='/home/ahsanm1/umair_wlab/data/NanoVar_data/pileups/training/chr%d/chr%d-' %(chrom,chrom)
                _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos',cpu=cpu)
                _,nx_train,ny_train,ntrain_allele,ntrain_ref= get_data(f_path+'neg',cpu=cpu)
                n_start=-len(x_train)
                n_end=0
                for i in range(10):

                    n_start+=len(x_train)
                    n_end=n_start+len(nx_train)
                    if n_end>len(x_train):
                        n_start=0
                        n_end=n_start+len(nx_train)
                    batch_nx_train,batch_ny_train,batch_ntrain_allele, batch_ntrain_ref = \
                    nx_train[n_start:n_end,:,:,:],ny_train[n_start:n_end,:],\
                    ntrain_allele[n_start:n_end,:], ntrain_ref[n_start:n_end,:]
                    
                    for batch in range(len(x_train)//batch_size):
                        batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))],batch_nx_train[batch*batch_size:min((batch+1)*batch_size,len(batch_nx_train))]])
                        batch_y = np.vstack([y_train[batch*batch_size:min((batch+1)*batch_size,len(y_train))],batch_ny_train[batch*batch_size:min((batch+1)*batch_size,len(batch_ny_train))]])    
                        batch_ref = np.vstack([train_ref[batch*batch_size:min((batch+1)*batch_size,len(train_ref))],batch_ntrain_ref[batch*batch_size:min((batch+1)*batch_size,len(batch_ntrain_ref))]])
                        batch_allele = np.vstack([train_allele[batch*batch_size:min((batch+1)*batch_size,len(train_allele))],batch_ntrain_allele[batch*batch_size:min((batch+1)*batch_size,len(batch_ntrain_allele))]])
                        # Run optimization op (backprop).
                            # Calculate batch loss and accuracy
                        opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele,keep:0.5})

                _,x_train,y_train,train_allele,train_ref= None,None,None,None,None
                _,nx_train,ny_train,ntrain_allele,ntrain_ref= None,None,None,None,None
                
            if k%5==4:
                    saver.save(sess, save_path=params['model'],global_step=save_num)
                    elapsed=time.time()-t
                    print ('Time Taken for Iteration %d-th: %.2f seconds\n'\
                           %((k+1)*10 ,elapsed), flush=True)
                    t=time.time()
                    save_num+=1
            
                
            
            

        #saver.save(sess, save_path=params['model'],global_step=save_num)
    

    
def genotype_caller_deprecated(params,input_type='path',data=None):
    tf.reset_default_graph()
    n_input=params['dims']
    if input_type=='path':
        _,x_train,y_train,train_allele,train_ref= get_data(params['train_path']+'pos')
        _,nx_train,ny_train,ntrain_allele,ntrain_ref= get_data(params['train_path']+'neg')
        print('train data received')
        _,px_test,py_test,ptest_allele,ptest_ref= get_data(params['test_path']+'pos')
        _,nx_test,ny_test,ntest_allele,ntest_ref= get_data(params['test_path']+'neg')
        x_test,y_test,test_allele,test_ref=np.vstack([px_test,nx_test]),np.vstack([py_test,ny_test]),\
        np.vstack([ptest_allele,ntest_allele]),np.vstack([ptest_ref,ntest_ref])
        print('test data received')

    else:
        _,x_train,y_train,train_allele,train_ref=data['train_data'][0]
        _,nx_train,ny_train,ntrain_allele,ntrain_ref= data['train_data'][1]
        _,px_test,py_test,ptest_allele,ptest_ref= data['test_data'][0]
        _,nx_test,ny_test,ntest_allele,ntest_ref= data['test_data'][1]
        x_test,y_test,test_allele,test_ref=np.vstack([px_test,nx_test]),np.vstack([py_test,ny_test]),\
        np.vstack([ptest_allele,ntest_allele]),np.vstack([ptest_ref,ntest_ref])

    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    n_input=[i for i in x_train.shape[1:]]

    weights,biases,t1,t2=get_tensors(n_input,learning_rate)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,rate)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2

    init = tf.global_variables_initializer()
    saver = tf.train.Saver()
    
    '''sample=np.random.choice(len(x_test), 300, replace=False)
    bx_test=x_test[sample]
    by_test=y_test[sample]
    btest_ref=test_ref[sample]
    btest_allele=test_allele[sample]'''
    
    n_size=1
    with tf.Session() as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        stats=[]
        print('starting training')
        n_start=-len(x_train)
        n_end=0
        for i in range(training_iters):
            n_start+=len(x_train)
            n_end+=len(x_train)
            if n_end>len(nx_train):
                n_start=0
                n_end=n_start+len(x_train)
            batch_nx_train,batch_ny_train,batch_ntrain_allele,batch_ntrain_ref=nx_train[n_start:n_end,:,:,:],\
                ny_train[n_start:n_end,:],ntrain_allele[n_start:n_end,:],ntrain_ref[n_start:n_end,:]
            
            for batch in range(len(x_train)//batch_size):
                
                batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))],batch_nx_train])
                batch_y = np.vstack([y_train[batch*batch_size:min((batch+1)*batch_size,len(y_train))],batch_ny_train])    
                batch_ref = np.vstack([train_ref[batch*batch_size:min((batch+1)*batch_size,len(train_ref))],batch_ntrain_ref])
                batch_allele = np.vstack([train_allele[batch*batch_size:min((batch+1)*batch_size,len(train_allele))],batch_ntrain_allele])
                # Run optimization op (backprop).
                    # Calculate batch loss and accuracy
                opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele,rate:0.5})
                loss, acc_gt,acc_allele = sess.run([cost, accuracy_gt,accuracy_allele], feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele,rate:0.5})
            print("Iter " + str(i) + ", Loss= " + \
                          "{:.6f}".format(loss) + ",\nTraining GT Accuracy= " + \
                          "{:.5f}".format(acc_gt)+",\nTraining Allele Accuracy= " + \
                          "{:.5f}".format(acc_allele))
            print("Optimization Finished!\n")

            # Calculate accuracy for all 10000 mnist test images
            sample=np.random.choice(len(x_test), 300, replace=False)
            bx_test=x_test[sample]
            by_test=y_test[sample]
            btest_ref=test_ref[sample]
            btest_allele=test_allele[sample]
            test_acc_gt,test_acc_allele,valid_loss = sess.run([accuracy_gt,accuracy_allele,cost],\
                                           feed_dict={x: bx_test,y : by_test,ref:btest_ref,allele:btest_allele,rate:0.0})
                # repeat steps 4 for the histogram summary            
            stats.append([loss,acc_gt,acc_allele,valid_loss,test_acc_gt,test_acc_allele])
            print("Testing Loss= "+"{:.5f}".format(valid_loss)+"\nTesting GT Accuracy= "+"{:.5f}".format(test_acc_gt)+",\nTesting Allele Accuracy= "+"{:.5f}".format(test_acc_allele))
            print('--------------------------------------------------\n')
        saver.save(sess, save_path=params['model'])
    stats=np.array(stats)
    np.save('stats_full',stats)
    return stats


def test_model(params,data=None,tp=None):
    model_path,test_path,n_input,chrom,vcf_path= params['model'], params['test_path'],params['dims'],params['chrom'],params['vcf_path']
    cpu=params['cpu']
    tf.reset_default_graph()
    if data:
        pos,x_test,test_ref=data
    else:
        
        pos,x_test,_,_,test_ref=get_data(params['test_path'],cpu=cpu,mode='test')
        
    print('test data received',flush=True)

    n_input=[i for i in x_test.shape[1:]]
    
    weights,biases,t1,t2=get_tensors(n_input,1)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2
    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    
    rev_mapping={0:'A',1:'G',2:'T',3:'C'}
    gt_map={0:1,1:0}
    
    batch_size=1000
    total=[]
    with open(vcf_path,'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        

        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus Genotype across all datasets with called genotype">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
        ttt=[]
        for batch in range(len(x_test)//batch_size+1):
                    batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                    batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]
                    
                    batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]
                    
                    # Run optimization op (backprop).
                        # Calculate batch loss and accuracy


                    fc_layer_batch,score_batch,gt_like,all_like = sess.run([fc_layer,pred,gt_likelihood,allele_likelihood],\
                                               feed_dict={x: batch_x,ref:batch_ref,keep:1.0})
                    
                    ref_like=np.max(all_like*batch_ref,axis=1)
                    qual=-10*np.log10(np.abs((gt_like[:,0]-1e-9)))-10*np.log10(np.abs((ref_like-1e-9)))
                    all_pred,gt_pred=np.argmax(fc_layer_batch,axis=1),np.argmax(score_batch,axis=1)
                    
                        
                    mat=np.hstack([batch_pos,np.argmax(batch_ref,axis=1)[:,np.newaxis],\
                               all_pred[:,np.newaxis],gt_pred[:,np.newaxis],qual[:,np.newaxis]])
                    total.append(mat)
                    mat=mat[(mat[:,1]!=mat[:,2])]

                    for j in range(len(mat)):

                        v=mat[j]

                        s='%s\t%d\t.\t%s\t%s\t%d\t%s\t.\tGT\t1/%d\n' %\
                        (chrom,v[0],rev_mapping[v[1]],rev_mapping[v[2]],v[4],'PASS',gt_map[v[3]])

                        f.write(s)
                            
                    
                              
    stats=np.array(total)
    np.save(vcf_path+'_stats',stats)
    return np.vstack(total)
                    
        
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
    
    args = parser.parse_args()
    input_dims=[int(x) for x in args.dimensions.split(':')]
    t=time.time()
    
    
    if args.mode=='train':
        in_dict={'cpu':args.cpu,'rate':args.rate,'iters':args.iterations,'size':args.size,'dims':input_dims,\
                 'train_path':args.train,'test_path':args.test,'model':args.model}
        genotype_caller_skinny(in_dict)
    
    else:
        in_dict={'cpu':args.cpu,'dims':input_dims,'test_path':args.test,'model':args.model,'chrom':args.chrom,'vcf_path':args.vcf}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
