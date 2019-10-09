import time,os,copy,argparse,subprocess,psutil
import pandas as pd
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect import *
from utils import *

config = tf.ConfigProto(device_count={"CPU": 64})
#config.gpu_options.allow_growth = True


def genotype_caller_skinny(params,input_type='path',data=None,attempt=0,neg_part='neg.combined'):
    tf.reset_default_graph()
    
    false_mltplr=3
    
    cpu=params['cpu']
    n_input=params['dims']
    dims=n_input
    chrom_list=[1]#list(range(2,23)) #params['chrom'].split(':') 
    #chrom_list=list(range(int(chrom_list[0]),int(chrom_list[1])+1))
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    weights,biases,tensors=get_tensors(n_input,learning_rate)
    (x,GT_label,A_label, G_label, T_label, C_label,GT_score, A_score, G_score, T_score, C_score, accuracy_GT, accuracy_A,  accuracy_G,  accuracy_T,  accuracy_C, prediction_accuracy_GT, prediction_accuracy_A,  prediction_accuracy_G,  prediction_accuracy_T,  prediction_accuracy_C, prediction_GT, prediction_A,  prediction_G,  prediction_T,  prediction_C, accuracy, cost, optimizer,  cost_GT, cost_A, cost_G, cost_T, cost_C,A_ref,G_ref,T_ref,C_ref,prob_GT,prob_A,prob_G,prob_T,prob_C,keep)=tensors

    if params['val']:
        val_list=[]
        for v_path in params['test_path'].split(':'):
            vx_test, vy_test, vtest_allele, vtest_ref=get_data_20plus(params)
            vtest_allele=make_allele(vy_test,vtest_allele,vtest_ref)
            val_list.append((v_path,(vx_test, vy_test, vtest_allele, vtest_ref)))
        
    
    init = tf.global_variables_initializer()
    saver = tf.train.Saver(max_to_keep=100)
    
    rec_size=14+dims[0]*dims[1]*dims[2]*6
    
    
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
                
             
                    
                _,x_train,y_train,train_allele,train_ref= get_data(f_path+'pos',cpu=cpu, dims=n_input)

                _,nx_train,ny_train,ntrain_allele,ntrain_ref=get_data(f_path+neg_part, cpu=cpu, dims=n_input)

                train_allele=make_allele(y_train,train_allele,train_ref)
                ntrain_allele=make_allele(ny_train,ntrain_allele,ntrain_ref)
                
                chrom_data['chr%d'%chrom]=((x_train,y_train,train_allele,train_ref),(nx_train,ny_train,ntrain_allele,ntrain_ref))
                print('chromosome %d done | '%chrom,end='',flush=True)

        print('Finished reading pileups',flush=True)
        
        for k in range(iter_steps):
            
            for chrom in chrom_list:
                print('Training on chrom %d ' %(chrom),end='',flush=True)
                (x_train,y_train,train_allele,train_ref),(nx_train,ny_train,ntrain_allele,ntrain_ref)=chrom_data['chr%d'%chrom]


                n_start=-false_mltplr*len(x_train)
                n_end=0
                for i in range(iters):


                    n_start+=false_mltplr*len(x_train)
                    n_end=n_start+false_mltplr*len(x_train)

                    if n_end>len(nx_train):
                        batch_nx_train= np.vstack([nx_train[n_start:,:,:,:],nx_train[:n_end-len(nx_train),:,:,:]]) 

                        batch_ny_train= np.vstack([ny_train[n_start:,:],ny_train[:n_end-len(nx_train),:]])
                        batch_ntrain_allele=np.vstack([ntrain_allele[n_start:,:],ntrain_allele[:n_end-len(nx_train),:]])
                        batch_ntrain_ref = np.vstack([ntrain_ref[n_start:,:],ntrain_ref[:n_end-len(nx_train),:]])


                        n_start=n_end-len(nx_train)
                        n_end=n_start+false_mltplr*len(x_train)
                    else:    
                        batch_nx_train,batch_ny_train,batch_ntrain_allele, batch_ntrain_ref = \
                        nx_train[n_start:n_end,:,:,:],ny_train[n_start:n_end,:],\
                        ntrain_allele[n_start:n_end,:], ntrain_ref[n_start:n_end,:]


                    training_loss=0
                    total_train_data=0

                    total_false_size=int(false_mltplr*batch_size)

                    for batch in range(int(np.ceil(len(x_train)/batch_size))):
                        batch_x = np.vstack([x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))],\
                                  batch_nx_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx_train))]])

                        batch_y = np.vstack([y_train[batch*batch_size:min((batch+1)*batch_size, len(y_train))],\
                                  batch_ny_train[total_false_size*batch :min(total_false_size*(batch+1),\
                                  len(batch_ny_train))]])

                        batch_allele = np.vstack([train_allele[batch*batch_size :min((batch+1)*batch_size,\
                                                   len(train_allele))], batch_ntrain_allele[total_false_size*batch : \
                                                   min(total_false_size*(batch+1), len(batch_ntrain_allele))]])

                        batch_ref = np.vstack([train_ref[batch*batch_size :min((batch+1)*batch_size,\
                                                len(train_ref))], batch_ntrain_ref[ total_false_size*batch:\
                                                min(total_false_size*(batch+1), len(batch_ntrain_ref))]])

                        opt,loss = sess.run([optimizer,cost], feed_dict={x: batch_x, GT_label:batch_y, A_label:np.eye(2)[batch_allele[:,0]], G_label:np.eye(2)[batch_allele[:,1]], T_label:np.eye(2)[batch_allele[:,2]], C_label:np.eye(2)[batch_allele[:,3]] , A_ref:batch_ref[:,0][:,np.newaxis], G_ref:batch_ref[:,1][:,np.newaxis], T_ref:batch_ref[:,2][:,np.newaxis], C_ref:batch_ref[:,3][:,np.newaxis], keep:0.5})
                        training_loss+=loss*len(batch_x)
                        total_train_data+=len(batch_x)


                    training_loss=training_loss/total_train_data

                    print('.',end='',flush=True)

                if params['val'] and (k<3 or chrom==22):

                    for val in val_list:
                        print('\n')
                        print(30*'-')
                        print(val[0])
                        vx_test, vy_test, vtest_allele, vtest_ref=val[1]


                        v_loss,v_acc,A_acc,G_acc,T_acc,C_acc,GT_acc,v_loss,v_acc=0,0,0,0,0,0,0,0,0

                        v_score,v_pred=[],[]

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

                        print('train loss= %.4f' %loss)
                        print('valid loss= %.4f\n' %(v_loss/len(vx_test)))
                        print('valid accuracy= %.4f' %(v_acc/len(vx_test)))
                        print('GT_acc=%.4f, A_acc=%.4f, G_acc=%.4f, T_acc=%.4f, C_acc=%.4f' %(GT_acc/len(vx_test), A_acc/len(vx_test), G_acc/len(vx_test), T_acc/len(vx_test), C_acc/len(vx_test)))
                        #print('Validation Precision= %.4f      Validation Recall= %.4f' %(true_positive/total_positive, true_positive/total_true))
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
    (x,GT_label,A_label, G_label, T_label, C_label,GT_score, A_score, G_score, T_score, C_score, accuracy_GT, accuracy_A,  accuracy_G,  accuracy_T,  accuracy_C, prediction_accuracy_GT, prediction_accuracy_A,  prediction_accuracy_G,  prediction_accuracy_T,  prediction_accuracy_C, prediction_GT, prediction_A,  prediction_G,  prediction_T,  prediction_C, accuracy, cost, optimizer,  cost_GT, cost_A, cost_G, cost_T, cost_C,A_ref,G_ref,T_ref,C_ref,prob_GT,prob_A,prob_G,prob_T,prob_C,keep)=tensors

    
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
    
    neg_file=open(vcf_path+'.neg'+suffix,'w')
    vcf_path=vcf_path+suffix
    
    with open(vcf_path,'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        

        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GP,Number=1,Type=Integer,Description="Genotype Probability">\n')
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
        total_prob=[]
        total_gt_prob=[]
        total_ref=[]
        for i in range(len(tmp_sz)-1):
            pos,x_test,_,_,test_ref= get_data(f_path,a=tmp_sz[i], b=tmp_sz[i+1],dims=n_input,cpu=cpu,mode='test')
            for batch in range(int(np.ceil(len(x_test)/batch_size))):
                        batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                        batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]

                        batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]

                        batch_prob_GT,batch_prob_A,batch_prob_G,batch_prob_C,batch_prob_T= sess.run([prob_GT,prob_A,prob_G,prob_C,prob_T],\
                                                   feed_dict={x: batch_x, A_ref:batch_ref[:,0][:,np.newaxis], G_ref:batch_ref[:,1][:,np.newaxis], T_ref:batch_ref[:,2][:,np.newaxis], C_ref:batch_ref[:,3][:,np.newaxis],keep:1.0})
                        
                        batch_pred_GT=np.argmax(batch_prob_GT,axis=1)
                        
                        batch_probs=np.hstack([batch_prob_A[:,1][:,np.newaxis], batch_prob_G[:,1][:,np.newaxis], batch_prob_T[:,1][:,np.newaxis], batch_prob_C[:,1][:,np.newaxis]])
                        
                        
                        total_prob.append(batch_probs)
                        
                        total_gt_prob.append(batch_prob_GT)
                        
                        batch_pred=np.argsort(batch_probs,axis=1)
                        
                        batch_ref=np.argmax(batch_ref,1)
                        total_ref.append(batch_ref[:,np.newaxis])
                        
                        batch_pred_GT=np.sum(batch_probs>=0.5,axis=1)
                        
                        for j in range(len(batch_pred_GT)):
                            
                            if batch_pred_GT[j]>=2: # if het
                                    pred1,pred2=batch_pred[j,-1],batch_pred[j,-2]
                                    if pred1==batch_ref[j] and batch_probs[j,pred2]>=0.99:
                                                s='%s\t%d\t.\t%s\t%s\t%d\t%s\t.\tGT:GP\t%s:%d\n' %(chrom, batch_pos[j], rev_mapping[batch_ref[j]], rev_mapping[pred2], 0,'PASS','0/1', int(min(99,-10*np.log10(1e-10+ 1-batch_probs[j,pred2]))))
                                                f.write(s)
                                        

                                    
                                    elif pred2==batch_ref[j] and batch_probs[j,pred2]>=0.99:
                                        s='%s\t%d\t.\t%s\t%s\t%d\t%s\t.\tGT:GP\t%s:%d\n' %(chrom,batch_pos[j], rev_mapping[batch_ref[j]], rev_mapping[pred1], 0,'PASS','1/0', int(min(99,-10*np.log10(1e-10+1-batch_probs[j,pred2]))))

                                        f.write(s)
                                        
                                    elif pred2!=batch_ref[j] and pred1!=batch_ref[j] and batch_probs[j,pred2]>=0.99:
                                        s='%s\t%d\t.\t%s\t%s,%s\t%d\t%s\t.\tGT:GP\t%s:%d\n' %\
                            (chrom,batch_pos[j],rev_mapping[batch_ref[j]],rev_mapping[pred1],rev_mapping[pred2],0,'PASS','1/2', int(min(99,-10*np.log10(1e-10+1-batch_probs[j,pred2]))))

                                        f.write(s)
                                
                            elif batch_pred_GT[j]==1 and batch_ref[j]!=batch_pred[j,-1] and batch_probs[j,batch_pred[j,-1]]>=0.99:
                                pred1=batch_pred[j,-1]
                                s='%s\t%d\t.\t%s\t%s\t%d\t%s\t.\tGT:GP\t%s:%d\n' %(chrom, batch_pos[j], rev_mapping[batch_ref[j]], rev_mapping[pred1], 0, 'PASS', '1/1', int(min(99,-10*np.log10(1e-10+1-batch_probs[j,pred1]))))
                                f.write(s)
                            else:
                                neg_file.write('%d,%.4f,%.4f,%.4f,%.4f,%.4f\n' %(batch_pos[j], batch_prob_GT[j,0], batch_probs[j,0], batch_probs[j,1], batch_probs[j,2], batch_probs[j,3]))
                                                                                                   
                                                                                                                   
                       
    neg_file.close()
    total_prob=np.array(total_prob)
    total_gt_prob=np.array(total_gt_prob)
    total_ref=np.array(total_ref)                                     
    if prob_save:
        np.savez('prob_file'+suffix,total_prob=total_prob,total_gt_prob=total_gt_prob,total_ref=total_ref)
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

def get_data_20plus(params):
    cpu=params['cpu']
    n_input=params['dims']
    f_path=params['test_path']
    _,vpx_train,vpy_train,vptrain_allele,vptrain_ref= get_data(f_path+'pos',cpu=cpu,dims=n_input)
    _,vnx_test,vny_test,vntest_allele,vntest_ref=get_data(f_path+'neg.15plus',cpu=cpu,dims=n_input)
    vx_test,vy_test,vtest_allele,vtest_ref =np.vstack([vpx_train,vnx_test]), np.vstack([vpy_train,vny_test]), np.vstack([vptrain_allele,vntest_allele]), np.vstack([vptrain_ref,vntest_ref])
    #v_truth=np.vstack((np.ones(len(vpx_train))[:,np.newaxis].astype(int), np.zeros(len(vnx_test))[:,np.newaxis].astype(int)))
    return (vx_test,vy_test,vtest_allele,vtest_ref)
    #return (vpx_train,vpy_train,vptrain_allele,vptrain_ref),(vnx_test,vny_test,vntest_allele,vntest_ref)

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
    
    args = parser.parse_args()
    input_dims=[int(x) for x in args.dimensions.split(':')]
    t=time.time()
    
    if args.mode=='train':
        in_dict={'cpu':args.cpu,'rate':args.rate, 'iters':args.iterations, 'size':args.size,'dims':input_dims,'chrom':args.chrom,\
                 'train_path':args.train, 'test_path':args.test, 'model':args.model, 'val':args.validation,'retrain':args.retrain,\
                'window':args.window,'ratio':args.ratio,'rt_path':args.rt_path}
        genotype_caller_skinny(in_dict,neg_part=args.neg_part)
    
    else:
        in_dict={'cpu':args.cpu,'dims':input_dims,'test_path':args.test,'model':args.model, 'chrom':args.chrom, 'vcf_path':args.vcf, 'workdir':args.workdir}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
