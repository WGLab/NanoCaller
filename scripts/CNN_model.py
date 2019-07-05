import time,os,copy,argparse,utils,subprocess,psutil
import pandas as pd
import numpy as np
import multiprocessing as mp
import generate_pileups as gcp
import tensorflow as tf


def get_data(fname,a=None, b=None,dims=(32,33,5), cpu=4,mode='train'):
    t=time.time()
    l=os.stat(fname).st_size
    
    if mode=='train':
        rec_size=14+dims[0]*dims[1]*7
        if a!=None and b!=None:
            my_array=[(fname,x,mode,dims) for x in range(a,b,1000*rec_size)]
        else:
            my_array=[(fname,x,mode,dims) for x in range(0,l,1000*rec_size)]
    else:
        rec_size=11+dims[0]*dims[1]*7
        my_array=[(fname,x,mode,dims) for x in range(a,b,1000*rec_size)]

    cpu=min(cpu,len(my_array))
    pool = mp.Pool(processes=cpu)
    results = pool.map(utils.read_pileups_from_file, my_array)
    pool.close()  
    pool.join() 
    
    
    pos=np.vstack([res[0][:,np.newaxis] for res in results])
    mat=np.vstack([res[1] for res in results])
    ref=np.vstack([res[2] for res in results])
    allele,gt=None,None
    if mode=='train':
        allele=np.vstack([res[3] for res in results])
        gt=np.vstack([res[4] for res in results])
    
    elapsed=time.time()-t
    
    print('I/O Time Elapsed: %.2f seconds' %elapsed, flush = True)
 
    
    return pos,mat,gt,allele,ref
        
def conv2d(x, W, b, strides=1,pad_type='VALID'):
    x = tf.nn.conv2d(x, W, strides=[1, strides, strides, 1], padding=pad_type)
    x = tf.nn.bias_add(x, b)
    return tf.nn.selu(x)




def conv_net(x,ref, weights, biases,keep):
    conv1_1 = conv2d(x, weights['wc1_1'], biases['bc1_1'],1,pad_type='SAME')
    conv1_2 = conv2d(x, weights['wc1_2'], biases['bc1_2'],1,pad_type='SAME')
    conv1_3 = conv2d(x, weights['wc1_3'], biases['bc1_3'],1,pad_type='SAME')
    
    merge_conv1=tf.concat([conv1_1, conv1_2,conv1_3],3)
    
    pooled=tf.nn.max_pool(merge_conv1,ksize=[1,4,1,1],strides=[1,1,1,1],padding='SAME')
    
    conv2 = conv2d(merge_conv1, weights['wc2'], biases['bc2'],2)
    
    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],2)

    flat_nn = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    
    drop_out_flat =flat_nn# tf.nn.dropout(flat_nn, keep)
    
    fc1 = tf.add(tf.matmul(drop_out_flat, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.selu(fc1)
    
    drop_out_fc1=tf.nn.dropout(fc1, keep)
    
    
    
    fc2 = tf.add(tf.matmul(drop_out_fc1, weights['wd-gtp']), biases['bd-gtp'])
    fc2 = tf.nn.selu(fc2)
    out_gt= tf.add(tf.matmul(fc2, weights['out-gtp']), biases['out-gtp'])
    
    
    fa = tf.add(tf.matmul(tf.concat([out_gt,drop_out_fc1,ref],1), weights['w-tr']), biases['w-tr'])
    fa = tf.nn.selu(fa)
    out_allele = tf.add(tf.matmul(fa, weights['out-all']), biases['out-all'])

    return out_gt,out_allele

def get_tensors(n_input,learning_rate=0):
    h_in,w_in,depth=n_input
    
    #h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    #w_in=int(np.ceil(float(w_in-4+1)/float(1)))
    
    for i in range(2):
        h_in=int(np.ceil(float(h_in-3+1)/float(2)))
        w_in=int(np.ceil(float(w_in-3+1)/float(2)))
       
    
    weights = {
        'wc1_1': tf.get_variable('W0_1', shape=(1,4,depth,16), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc1_2': tf.get_variable('W0_2', shape=(32,1,depth,16), initializer=tf.contrib.layers.xavier_initializer()),
        'wc1_3': tf.get_variable('W0_3', shape=(3,3,depth,16), initializer=tf.contrib.layers.xavier_initializer()),
        'wc2': tf.get_variable('W1', shape=(3,3,48,64), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(3,3,64,128), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(h_in*w_in*128,64), initializer=tf.contrib.layers.xavier_initializer()), 
        'w-tr': tf.get_variable('W4', shape=(70,32), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('W5', shape=(32,4), initializer=tf.contrib.layers.xavier_initializer()),
        'wd-gtp': tf.get_variable('W6', shape=(64,32), initializer=tf.contrib.layers.xavier_initializer()), 
        'out-gtp': tf.get_variable('W7', shape=(32,2), initializer=tf.contrib.layers.xavier_initializer()), 
    }
    
    biases = {
        'bc1_1': tf.get_variable('B0_1', shape=(16), initializer=tf.contrib.layers.xavier_initializer()), 
        'bc1_2': tf.get_variable('B0_2', shape=(16), initializer=tf.contrib.layers.xavier_initializer()),
        'bc1_3': tf.get_variable('B0_3', shape=(16), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(64), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(64), initializer=tf.contrib.layers.xavier_initializer()),
        'w-tr': tf.get_variable('B4', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('B5', shape=(4), initializer=tf.contrib.layers.xavier_initializer()),
        'bd-gtp': tf.get_variable('B6', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'out-gtp': tf.get_variable('B7', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
    }
    
    '''weights = {
        'wc1': tf.get_variable('W0', shape=(conv_krnl_shapes[0][0], conv_krnl_shapes[0][1] ,depth,32), initializer=tf.contrib.layers.variance_scaling_initializer()), 
        
        'wc2': tf.get_variable('W1', shape=(conv_krnl_shapes[1][0], conv_krnl_shapes[1][1],32,64), initializer=tf.contrib.layers.variance_scaling_initializer()), 
        
        'wc3': tf.get_variable('W2', shape=(conv_krnl_shapes[2][0], conv_krnl_shapes[2][1] ,64,128), initializer=tf.contrib.layers.variance_scaling_initializer()), 
        
        'wd1': tf.get_variable('W3', shape=(conv_shapes[-1][0]*conv_shapes[-1][1]*128, 64), initializer= tf.contrib.layers.variance_scaling_initializer()), 
        
        'w-tr': tf.get_variable('W4', shape=(100,32), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'out-all': tf.get_variable('W5', shape=(32,4), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'wd-gtp': tf.get_variable('W6', shape=(64,32), initializer=tf.contrib.layers.variance_scaling_initializer()), 
        'out-gtp': tf.get_variable('W7', shape=(32,2), initializer=tf.contrib.layers.variance_scaling_initializer()), 
    }
    
    biases = {
        'bc1': tf.get_variable('B0', shape=(32), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'bc2': tf.get_variable('B1', shape=(64), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'bc3': tf.get_variable('B2', shape=(128), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'bd1': tf.get_variable('B3', shape=(64), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'w-tr': tf.get_variable('B4', shape=(32), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'out-all': tf.get_variable('B5', shape=(4), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'bd-gtp': tf.get_variable('B6', shape=(32), initializer=tf.contrib.layers.variance_scaling_initializer()),
        'out-gtp': tf.get_variable('B7', shape=(2), initializer=tf.contrib.layers.variance_scaling_initializer()),
    }'''
    
    x,y = tf.placeholder("float", [None]+n_input), tf.placeholder("float", [None, 2])
    allele,ref=tf.placeholder("float", [None, 4]),tf.placeholder("float", [None, 4])
    keep = tf.placeholder(tf.float32)
    
    pred,fc_layer = conv_net(x,ref,weights, biases, keep)
    
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
        f_path=params['train_path']
        _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos',cpu=cpu,dims=n_input)

        negative_variants=[get_data(f_path+'neg_%d' %freq,cpu=cpu,dims=n_input) for freq in [0,5,10,15,20,25]]

        nx_train=np.vstack([tmp[1] for tmp in negative_variants])
        ny_train=np.vstack([tmp[2] for tmp in negative_variants])
        ntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
        ntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

        perm=np.random.permutation(len(nx_train))

        np.take(nx_train,perm,axis=0,out=nx_train)
        np.take(ny_train,perm,axis=0,out=ny_train)
        np.take(ntrain_allele,perm,axis=0,out=ntrain_allele)
        np.take(ntrain_ref,perm,axis=0,out=ntrain_ref)

        f_path=params['test_path']
        _,vpx_train,vpy_train,vptrain_allele,vptrain_ref= get_data(f_path+'pos',cpu=cpu,dims=n_input)

        negative_variants=[get_data(f_path+'neg_%d' %freq,cpu=cpu,dims=n_input) for freq in [0,5,10,15,20,25]]

        vnx_train=np.vstack([tmp[1] for tmp in negative_variants])
        vny_train=np.vstack([tmp[2] for tmp in negative_variants])
        vntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
        vntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

        perm=np.random.permutation(len(vnx_train))

        np.take(vnx_train,perm,axis=0,out=vnx_train)
        np.take(vny_train,perm,axis=0,out=vny_train)
        np.take(vntrain_allele,perm,axis=0,out=vntrain_allele)
        np.take(vntrain_ref,perm,axis=0,out=vntrain_ref)
        
        test_size=len(vpx_train)//10
        vx_test,vy_test,vtest_allele,vtest_ref=np.vstack([vpx_train,vnx_train[:test_size]]),\
        np.vstack([vpy_train,vny_train[:test_size]]),np.vstack([vptrain_allele,vntrain_allele[:test_size]]),\
        np.vstack([vptrain_ref,vntrain_ref[:test_size]])
        print('train data received')
    else:
        x_train,y_train,train_allele,train_ref=data['train_data'][0]
        nx_train,ny_train,ntrain_allele,ntrain_ref= data['train_data'][1]
        
        vx_test,vy_test,vtest_allele,vtest_ref=data['test_data']
        
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
            print('iteration number: %d' %i,flush=True)
            false_mltplr=2
            n_start+=false_mltplr*len(x_train)
            n_end=n_start+false_mltplr*len(nx_train)
            if n_end>len(x_train):
                n_start=0
                n_end=n_start+false_mltplr*len(nx_train)
            batch_nx_train,batch_ny_train,batch_ntrain_allele, batch_ntrain_ref = \
            nx_train[n_start:n_end,:,:,:],ny_train[n_start:n_end,:],\
            ntrain_allele[n_start:n_end,:], ntrain_ref[n_start:n_end,:]

            for batch in range(len(x_train)//batch_size):
                batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))],\
                          batch_nx_train[ false_mltplr*batch*batch_size : min(false_mltplr*(batch+1)*batch_size,\
                          len(batch_nx_train))]])

                batch_y = np.vstack([y_train[batch*batch_size:min((batch+1)*batch_size, len(y_train))],\
                          batch_ny_train[false_mltplr*batch*batch_size :min(false_mltplr*(batch+1)*batch_size,\
                          len(batch_ny_train))]])    
                batch_ref = np.vstack([train_ref[batch*batch_size :min((batch+1)*batch_size,\
                            len(train_ref))], batch_ntrain_ref[ false_mltplr*batch*batch_size:\
                            min(false_mltplr*(batch+1)*batch_size, len(batch_ntrain_ref))]])

                batch_allele = np.vstack([train_allele[batch*batch_size :min((batch+1)*batch_size,\
                               len(train_allele))], batch_ntrain_allele[false_mltplr*batch*batch_size : \
                               min(false_mltplr*(batch+1)*batch_size, len(batch_ntrain_allele))]])
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
            print('train loss= %.4f   train loss diff= %.4f\n' %(loss,abs(loss_list[-1]-loss_list[-2])),flush=True)

            test_mat=[]
            for batch in range(len(vx_test)//(batch_size)):
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
    chrom_list=[int(x) for x in params['chrom'].split(':')]
    
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
    with tf.Session()  as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        stats=[]
        v_stats=[]
        print('starting training',flush=True)
        count=0
        
        save_num=1
        t=time.time()
        
        iter_ratio=50
        
        iter_steps=max(training_iters//iter_ratio,1)
        iters=min(iter_ratio,training_iters)
        
        for k in range(iter_steps):
            
            for chrom in chrom_list:
                if chrom<10:
                    chnk=10
                elif chrom<16:
                    chnk=8
                else:
                    chnk=4
                tot_list={}
                f_path=os.path.join(params['train_path'],'chr%d/chr%d_pileups_' %(chrom,chrom))
                statinfo = os.stat(f_path+'pos')
                sz=statinfo.st_size
                tmp_sz=list(range(0,sz,7406000*(sz//(chnk*7406000))))
                tmp_sz=tmp_sz[:chnk]
                tmp_sz=tmp_sz+[sz] if tmp_sz[-1]!=sz else tmp_sz
                tot_list['pos']=tmp_sz


                for i in [0,5,10,15,20,25]:
                    statinfo = os.stat(f_path+'neg_%d' %i)
                    sz=statinfo.st_size
                    tmp_sz=list(range(0,sz,7406000*(sz//(chnk*7406000))))
                    tmp_sz=tmp_sz[:chnk]
                    tmp_sz=tmp_sz+[sz] if tmp_sz[-1]!=sz else tmp_sz
                    tot_list[i]=tmp_sz


                for i in range(len(tot_list['pos'])-1):
                    _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos',a=tot_list['pos'][i], b=tot_list['pos'][i+1], cpu=cpu, dims=n_input)
                    negative_variants=[get_data(f_path+'neg_%d' %freq, a=tot_list[freq][i], b=tot_list[freq][i+1], cpu=cpu, dims=n_input) for freq in [0,5,10,15,20,25]]
                    
                    nx_train=np.vstack([tmp[1] for tmp in negative_variants])
                    ny_train=np.vstack([tmp[2] for tmp in negative_variants])
                    ntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
                    ntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

                    perm=np.random.permutation(len(nx_train))

                    np.take(nx_train,perm,axis=0,out=nx_train)
                    np.take(ny_train,perm,axis=0,out=ny_train)
                    np.take(ntrain_allele,perm,axis=0,out=ntrain_allele)
                    np.take(ntrain_ref,perm,axis=0,out=ntrain_ref)
                    false_mltplr=2
                    n_start=-false_mltplr*len(x_train)
                    n_end=0
                    for i in range(iters):

                        n_start+=false_mltplr*len(x_train)
                        n_end=n_start+false_mltplr*len(nx_train)
                        if n_end>len(x_train):
                            n_start=0
                            n_end=n_start+false_mltplr*len(nx_train)
                        batch_nx_train,batch_ny_train,batch_ntrain_allele, batch_ntrain_ref = \
                        nx_train[n_start:n_end,:,:,:],ny_train[n_start:n_end,:],\
                        ntrain_allele[n_start:n_end,:], ntrain_ref[n_start:n_end,:]

                        for batch in range(len(x_train)//batch_size):
                            batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))],\
                                      batch_nx_train[ false_mltplr*batch*batch_size : min(false_mltplr*(batch+1)*batch_size,\
                                      len(batch_nx_train))]])
                            
                            batch_y = np.vstack([y_train[batch*batch_size:min((batch+1)*batch_size, len(y_train))],\
                                      batch_ny_train[false_mltplr*batch*batch_size :min(false_mltplr*(batch+1)*batch_size,\
                                      len(batch_ny_train))]])    
                            batch_ref = np.vstack([train_ref[batch*batch_size :min((batch+1)*batch_size,\
                                        len(train_ref))], batch_ntrain_ref[ false_mltplr*batch*batch_size:\
                                        min(false_mltplr*(batch+1)*batch_size, len(batch_ntrain_ref))]])
                            
                            batch_allele = np.vstack([train_allele[batch*batch_size :min((batch+1)*batch_size,\
                                           len(train_allele))], batch_ntrain_allele[false_mltplr*batch*batch_size : \
                                           min(false_mltplr*(batch+1)*batch_size, len(batch_ntrain_allele))]])
                            # Run optimization op (backprop).
                                # Calculate batch loss and accuracy
                            opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele,keep:0.5})
                        #utils.gpu_stats()
                    _,x_train,y_train,train_allele,train_ref= None,None,None,None,None
                    _,nx_train,ny_train,ntrain_allele,ntrain_ref= None,None,None,None,None
                    print('next chunk for chrom %d' %chrom,flush=True)
                    utils.usage_stats()
            saver.save(sess, save_path=params['model'],global_step=save_num)
            elapsed=time.time()-t
            
            print ('Time Taken for Iteration %d-%d: %.2f seconds\n'\
                   %((save_num-1)*iters,save_num*iters,elapsed), flush=True)
            
            utils.usage_stats()
            save_num+=1
            t=time.time()
        #saver.save(sess, save_path=params['model'],global_step=save_num)
    


def test_model(params,data=None,tp=None):
    model_path,test_path,n_input,chrom,vcf_path= params['model'], params['test_path'],params['dims'],params['chrom'],params['vcf_path']
    cpu=params['cpu']
    tf.reset_default_graph()
    '''if data:
        pos,x_test,test_ref=data
    else:
        
        pos,x_test,_,_,test_ref=get_data(params['test_path'],cpu=cpu,mode='test')'''

    
    
    
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
        
        chnk=15
        tot_list={}
        f_path=params['test_path']
        
        statinfo = os.stat(f_path)
        sz=statinfo.st_size
        tmp_sz=list(range(0,sz,7404000*(sz//(chnk*7404000))))
        tmp_sz=tmp_sz[:chnk]
        tmp_sz=tmp_sz+[sz] if tmp_sz[-1]!=sz else tmp_sz

        for i in range(len(tmp_sz)-1):
            pos,x_test,_,_,test_ref= get_data(f_path,a=tmp_sz[i], b=tmp_sz[i+1],cpu=cpu,mode='test')
            utils.usage_stats()
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
                    
    
def test_model_novcf(params,input_type='path',data=None):
    tf.reset_default_graph()
    cpu=params['cpu']
    n_input=params['dims']
    if input_type=='path':
        f_path=params['test_path']
        _,vpx_train,vpy_train,vptrain_allele,vptrain_ref= get_data(f_path+'pos',cpu=cpu)

        negative_variants=[get_data(f_path+'neg_%d' %freq,cpu=cpu) for freq in [0,5,10,15,20,25]]

        vnx_train=np.vstack([tmp[1] for tmp in negative_variants])
        vny_train=np.vstack([tmp[2] for tmp in negative_variants])
        vntrain_allele=np.vstack([tmp[3] for tmp in negative_variants])
        vntrain_ref=np.vstack([tmp[4] for tmp in negative_variants])

        perm=np.random.permutation(len(vnx_train))

        np.take(vnx_train,perm,axis=0,out=vnx_train)
        np.take(vny_train,perm,axis=0,out=vny_train)
        np.take(vntrain_allele,perm,axis=0,out=vntrain_allele)
        np.take(vntrain_ref,perm,axis=0,out=vntrain_ref)
        
        test_size=len(vpx_train)
        vx_test,vy_test,vtest_allele,vtest_ref=np.vstack([vpx_train,vnx_train[:test_size]]),\
        np.vstack([vpy_train,vny_train[:test_size]]),np.vstack([vptrain_allele,vntrain_allele[:test_size]]),\
        np.vstack([vptrain_ref,vntrain_ref[:test_size]])
        print('train data received')
    else:
        vx_test,vy_test,vtest_allele,vtest_ref=data['test_data']
        
    training_iters, learning_rate,model_path= params['iters'],\
    params['rate'], params['model']

    weights,biases,t1,t2=get_tensors(n_input,1)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2
    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    
    batch_size=1000
    total_stats={'num':0,'acc':0,'gt':0,'allele':0}
    for batch in range(len(vx_test)//(batch_size)):
        vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
        vbatch_y = vy_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))] 
        vbatch_ref = vtest_ref[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
        vbatch_allele = vtest_allele[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]

        fc_layer_batch,score_batch,v_loss = sess.run([fc_layer,pred,cost], feed_dict={x: vbatch_x,y: vbatch_y,ref:vbatch_ref,allele:vbatch_allele,keep:1.0})

        all_pred,gt_pred=np.argmax(fc_layer_batch,axis=1),np.argmax(score_batch,axis=1)

        test_mat=np.hstack([all_pred[:,np.newaxis],np.argmax(vbatch_allele,axis=1)[:,np.newaxis],\
                                   np.argmax(vbatch_ref,axis=1)[:,np.newaxis], gt_pred[:,np.newaxis],\
                                   np.argmax(vbatch_y,axis=1)[:,np.newaxis]])
        total_stats['num']+=len(test_mat)
        total_stats['acc']+=sum((test_mat[:,0]==test_mat[:,1])&(test_mat[:,4]==test_mat[:,3]))
        total_stats['gt']+=sum(test_mat[:,4]==test_mat[:,3])
        total_stats['allele']+=sum(test_mat[:,0]==test_mat[:,1])
        
    print('valid accuracy= %.4f' %(total_stats['acc']/total_stats['num']))
    print('valid GT accuracy= %.4f' %(total_stats['gt']/total_stats['num']))
    print('valid Allele accuracy= %.4f' %(total_stats['allele']/total_stats['num']))
   
    

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
        in_dict={'cpu':args.cpu,'rate':args.rate,'iters':args.iterations,'size':args.size,'dims':input_dims,'chrom':args.chrom,\
                 'train_path':args.train,'test_path':args.test,'model':args.model}
        genotype_caller_skinny(in_dict)
    
    else:
        in_dict={'cpu':args.cpu,'dims':input_dims,'test_path':args.test,'model':args.model,'chrom':args.chrom,'vcf_path':args.vcf}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
