import time,os,copy,argparse,subprocess,psutil
import pandas as pd
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect import *
from utils import *

config = tf.ConfigProto()
config.gpu_options.allow_growth = True

def genotype_caller(params,input_type='path',data=None):
    tf.reset_default_graph()
    cpu=params['cpu']
    n_input=params['dims']
    if input_type=='path':
        pos,neg=get_train_test(params,verbose=True)
        _,x_train,y_train,train_allele,train_ref= pos
        _,nx_train,ny_train,ntrain_allele,ntrain_ref=neg
        vx_test,vy_test,vtest_allele,vtest_ref=get_train_test(params,mode='test',verbose=True)
        
        sm_x_train,sm_y_train,sm_train_allele,sm_train_ref= x_train[len(x_train)//2: int(0.51*len(x_train))], y_train[len(x_train)//2: int(0.51*len(x_train))], train_allele[len(x_train)//2:int(0.51*len(x_train))], train_ref[len(x_train)//2:int(0.51*len(x_train))]
        sm_nx_train,sm_ny_train,sm_ntrain_allele,sm_ntrain_ref=nx_train[len(nx_train)//2: int(0.51*len(nx_train))], ny_train[len(nx_train)//2: int(0.51*len(nx_train))], ntrain_allele[len(nx_train)//2:int(0.51*len(nx_train))], ntrain_ref[len(nx_train)//2:int(0.51*len(nx_train))]
        
        sm_x_test=np.vstack([sm_x_train,sm_nx_train])
        sm_y_test=np.vstack([sm_y_train,sm_ny_train])
        sm_test_allele=np.vstack([sm_train_allele,sm_ntrain_allele])
        sm_test_ref=np.vstack([sm_train_ref,sm_ntrain_ref])
        
        
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
    (correct_prediction, correct_prediction_gt, correct_prediction_allele, accuracy, accuracy_gt, accuracy_allele, gt_likelihood, allele_likelihood)=t2
    
    
    
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
            test_stats={'num':0,'acc':0,'gt':0,'allele':0}
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
            tot_num,acc, acc_gt,acc_allele=0,0,0,0
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
            
            batch=(len(x_train)//batch_size)//2
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
            loss,acc, acc_gt,acc_allele,fc_layer_batch,score_batch = sess.run([cost, accuracy ,accuracy_gt, accuracy_allele, fc_layer,pred], feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele,keep:1.0})

            acc, acc_gt,acc_allele=acc/len(batch_x), acc_gt/len(batch_x),acc_allele/len(batch_x)
            v_error=1-acc
            loss_list.append(loss)
            loss_diff=max(loss_list[-p:])-min(loss_list[-p:])
            tp,fp,true=0,0,0
            
            for batch in range(len(vx_test)//(batch_size)):
                vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                vbatch_y = vy_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))] 
                vbatch_ref = vtest_ref[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                vbatch_allele = vtest_allele[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                
                
                fc_layer_batch,score_batch,v_loss,v_acc,v_gt_acc,v_all_acc,prediction = sess.run([fc_layer, pred, cost, accuracy, accuracy_gt, accuracy_allele, correct_prediction], feed_dict={x: vbatch_x,y: vbatch_y,ref:vbatch_ref, allele:vbatch_allele,keep:1.0})
                
                mat=np.hstack([prediction[:,np.newaxis], np.argmax(vbatch_y,axis=1)[:,np.newaxis],\
                                   np.argmax(vbatch_ref,axis=1)[:,np.newaxis], np.argmax(vbatch_allele,axis=1)[:,np.newaxis]])

                tmp=mat[mat[:,2]!=mat[:,3]]
                tp+=np.sum(tmp[:,0])
                true+=len(mat[mat[:,2]!=mat[:,3]])
                tmp=mat[mat[:,2]==mat[:,3]]
                fp+=(len(mat[mat[:,2]==mat[:,3]])-np.sum(tmp[:,0]))
                
                

                test_stats['num']+=len(vbatch_x)
                test_stats['acc']+=v_acc
                test_stats['gt']+=v_gt_acc
                test_stats['allele']+=v_all_acc
            print('train loss= %.4f   train loss diff= %.4f   valid loss= %.4f\n' %(loss,abs(loss_list[-1]-loss_list[-2]), v_loss), flush=True)
           
            
            print('train accuracy= %.4f           valid accuracy= %.4f' %(acc,test_stats['acc']/test_stats['num']))
            print('train GT accuracy= %.4f        valid GT accuracy= %.4f' %(acc_gt,test_stats['gt']/test_stats['num']))
            print('train Allele accuracy= %.4f    valid Allele accuracy= %.4f' %(acc_allele,test_stats['allele']/test_stats['num']))
            print('Validation Precision= %.4f      Validation Recall= %.4f' %(tp/(tp+fp),tp/true))
            print(100*'.')
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
    
    return


def genotype_caller_skinny(params,input_type='path',data=None,attempt=0):
    tf.reset_default_graph()
    cpu=params['cpu']
    n_input=params['dims']
    
    chrom_list=list(range(2,23))
    
    
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    weights,biases,t1,t2=get_tensors(n_input,learning_rate)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)=t1
    (correct_prediction, correct_prediction_gt, correct_prediction_allele, accuracy, accuracy_gt, accuracy_allele, gt_likelihood, allele_likelihood)=t2

    if params['val']:
        val_list=[]
        for v_path in params['test_path'].split(':'):
            val_list.append((v_path,get_train_test(params, mode='test', path=v_path)))
        
    
    init = tf.global_variables_initializer()
    saver = tf.train.Saver(max_to_keep=100)
    '''if params['retrain']:
        saver.restore(sess, model_path)
    restart=0'''
    rec_size=1000*(14+n_input[0]*n_input[1]*7)
    
    
    n_size=1
    with tf.Session(config=config)  as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        stats=[]
        v_stats=[]
        print('starting training',flush=True)
        count=0
        
        save_num=1
        t=time.time()
        
        iter_ratio=20
        
        iter_steps=max(training_iters//iter_ratio,1)
        iters=min(iter_ratio,training_iters)
        
        for k in range(iter_steps):
            
            for chrom in chrom_list:
                print('Training on chrom %d ' %(chrom),end='',flush=True)
               

                f_path=os.path.join(params['train_path'],'chr%d/chr%d.pileups.' %(chrom,chrom))
                
                _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos', cpu=params['cpu'], dims=params['dims'])

                _,nx_train,ny_train,ntrain_allele,ntrain_ref= get_data(f_path+'neg.combined', cpu=params['cpu'], dims=params['dims'])

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
                        batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size, len(x_train))],\
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

                print('.',end='',flush=True)

            if params['val']:

                for val in val_list:
                    print('\n')
                    print(30*'-')
                    print(val[0])
                    vx_test, vy_test, vtest_allele, vtest_ref=val[1]
                    test_stats={'num':0,'acc':0,'gt':0,'allele':0}
                    tp,true,fp=0,0,0
                    loss = sess.run(cost, feed_dict={x: batch_x,y: batch_y, ref:batch_ref, allele:batch_allele, keep:1})
                    for batch in range(len(vx_test)//(batch_size)):
                        vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                        vbatch_y = vy_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))] 
                        vbatch_ref = vtest_ref[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
                        vbatch_allele = vtest_allele[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]


                        fc_layer_batch,score_batch,v_loss,v_acc,v_gt_acc,v_all_acc,prediction = sess.run([fc_layer, pred, cost, accuracy, accuracy_gt, accuracy_allele, correct_prediction], feed_dict={x: vbatch_x,y: vbatch_y,ref:vbatch_ref, allele:vbatch_allele,keep:1.0})

                        mat=np.hstack([prediction[:,np.newaxis], np.argmax(vbatch_y,axis=1)[:,np.newaxis],\
                                   np.argmax(vbatch_ref,axis=1)[:,np.newaxis], np.argmax(vbatch_allele,axis=1)[:,np.newaxis]])
                        tmp=mat[mat[:,2]!=mat[:,3]]
                        tp+=np.sum(tmp[:,0])
                        true+=len(mat[mat[:,2]!=mat[:,3]])
                        tmp=mat[mat[:,2]==mat[:,3]]
                        fp+=(len(mat[mat[:,2]==mat[:,3]])-np.sum(tmp[:,0]))

                        test_stats['num']+=len(vbatch_x)
                        test_stats['acc']+=v_acc
                        test_stats['gt']+=v_gt_acc
                        test_stats['allele']+=v_all_acc


                    print('training loss= %.4f     valid loss= %.4f\n' %(loss, v_loss), flush=True)
                    print('valid accuracy= %.4f' %(test_stats['acc']/test_stats['num']), flush=True)
                    print('valid GT accuracy= %.4f' %(test_stats['gt']/test_stats['num']), flush=True)
                    print('valid Allele accuracy= %.4f' %(test_stats['allele']/test_stats['num']), flush=True)
                    print('validation Precision= %.4f     Validation Recall= %.4f' %(tp/(tp+fp),tp/true), flush=True)
                    print(30*'-')
                    print('\n')

            saver.save(sess, save_path=params['model'],global_step=save_num)
            elapsed=time.time()-t
            
            print ('Time Taken for Iteration %d-%d: %.2f seconds\n'\
                   %((save_num-1)*iters,save_num*iters,elapsed), flush=True)
            
            save_num+=1
            t=time.time()
            
        #saver.save(sess, save_path=params['model'],global_step=save_num)

def test_model(params):
    model_path,test_path,n_input,chrom,vcf_path= params['model'], params['test_path'],params['dims'],params['chrom'],params['vcf_path']
    cpu=params['cpu']
    tf.reset_default_graph()
    
    tr_dim=n_input[:]
    params['window']=None

    weights,biases,t1,t2=get_tensors(tr_dim,1)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)=t1
    (correct_prediction, correct_prediction_gt, correct_prediction_allele, accuracy, accuracy_gt, accuracy_allele, gt_likelihood, allele_likelihood)=t2
    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    
    rev_mapping={0:'A',1:'G',2:'T',3:'C'}
    gt_map={0:1,1:0}
    rec_size=1000*(12+n_input[0]*n_input[1]*7)
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
        
        chnk=1
        tot_list={}
        f_path=params['test_path']
        
        statinfo = os.stat(f_path)
        sz=statinfo.st_size
        tmp_sz=list(range(0,sz,rec_size*(sz//(chnk*rec_size))))
        tmp_sz=tmp_sz[:chnk]
        tmp_sz=tmp_sz+[sz] if tmp_sz[-1]!=sz else tmp_sz

        for i in range(len(tmp_sz)-1):
            pos,x_test,_,_,test_ref= get_data(f_path,a=tmp_sz[i], b=tmp_sz[i+1],dims=n_input,cpu=cpu,mode='test')
            for batch in range(len(x_test)//batch_size+1):
                        batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                        batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]

                        batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]

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
    return
                    
    
def test_model_novcf(params,input_type='path',data=None):
    tf.reset_default_graph()
    cpu=params['cpu']
    n_input=params['dims']
    if input_type=='path':
        vx_test,vy_test,vtest_allele,vtest_ref=get_train_test(params,mode='test')
    else:
        vx_test,vy_test,vtest_allele,vtest_ref=data['test_data']
        
    training_iters, learning_rate,model_path= params['iters'],\
    params['rate'], params['model']

    weights,biases,t1,t2=get_tensors(n_input,1)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)=t1
    (correct_prediction, correct_prediction_gt, correct_prediction_allele, accuracy, accuracy_gt, accuracy_allele, gt_likelihood, allele_likelihood)=t2
    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    
    batch_size=1000
    tp,fp,true=0,0,0
    test_stats={'num':0,'acc':0,'gt':0,'allele':0}
    for batch in range(len(vx_test)//(batch_size)):
        vbatch_x = vx_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
        vbatch_y = vy_test[batch*batch_size:min((batch+1)*batch_size,len(vx_test))] 
        vbatch_ref = vtest_ref[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]
        vbatch_allele = vtest_allele[batch*batch_size:min((batch+1)*batch_size,len(vx_test))]

        fc_layer_batch,score_batch,v_loss,v_acc,v_gt_acc,v_all_acc,prediction = sess.run([fc_layer, pred, cost, accuracy, accuracy_gt, accuracy_allele, correct_prediction], feed_dict={x: vbatch_x,y: vbatch_y,ref:vbatch_ref, allele:vbatch_allele,keep:1.0})

        mat=np.hstack([prediction[:,np.newaxis], np.argmax(vbatch_y,axis=1)[:,np.newaxis],\
                   np.argmax(vbatch_ref,axis=1)[:,np.newaxis], np.argmax(vbatch_allele,axis=1)[:,np.newaxis]])
        tmp=mat[mat[:,2]!=mat[:,3]]
        tp+=np.sum(tmp[:,0])
        true+=len(mat[mat[:,2]!=mat[:,3]])
        tmp=mat[mat[:,2]==mat[:,3]]
        fp+=(len(mat[mat[:,2]==mat[:,3]])-np.sum(tmp[:,0]))

        test_stats['num']+=len(vbatch_x)
        test_stats['acc']+=v_acc
        test_stats['gt']+=v_gt_acc
        test_stats['allele']+=v_all_acc

    print(100*'.')
    print('valid loss= %.4f\n' %( v_loss), flush=True)
    print('valid accuracy= %.4f' %(test_stats['acc']/test_stats['num']), flush=True)
    print('valid GT accuracy= %.4f' %(test_stats['gt']/test_stats['num']), flush=True)
    print('valid Allele accuracy= %.4f' %(test_stats['allele']/test_stats['num']), flush=True)
    print(' Validation Precision= %.4f     Validation Recall= %.4f' %(tp/(tp+fp),tp/true), flush=True)
    print(100*'.')
    print('\n')

   
    

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
    
    args = parser.parse_args()
    input_dims=[int(x) for x in args.dimensions.split(':')]
    t=time.time()
    
    if args.mode=='train':
        in_dict={'cpu':args.cpu,'rate':args.rate, 'iters':args.iterations, 'size':args.size,'dims':input_dims,'chrom':args.chrom,\
                 'train_path':args.train, 'test_path':args.test, 'model':args.model, 'val':args.validation,'retrain':args.retrain,\
                'window':args.window}
        genotype_caller_skinny(in_dict)
    
    else:
        in_dict={'cpu':args.cpu,'dims':input_dims,'test_path':args.test,'model':args.model,'chrom':args.chrom,'vcf_path':args.vcf}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
