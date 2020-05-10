import time,os,copy,argparse,subprocess,psutil
import pandas as pd
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect_indel_real_new_gt import *
from matplotlib import pyplot as plt
import git
from Bio import SeqIO
from Bio.Align.Applications import TCoffeeCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

config =  tf.compat.v1.ConfigProto(device_count={"CPU": 64})
#config.gpu_options.allow_growth = True

rev_gt_map={0:'hom-ref', 1:'hom-alt', 2:'het-ref', 3:'het-alt'}
rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

def pairwise(x,y):
    alignments = pairwise2.align.localms(x, y, 2, -1.0, -0.5, -0.1)

    # Use format_alignment method to format the alignments in the list
    #center=''

    #print(format_alignment(*alignments[0]))

    return alignments

def genotype_caller_skinny(params, input_type='path', data=None, neg_part= 'neg.combined'):
    tf.reset_default_graph()
    
    false_mltplr=1
    
    cpu=params['cpu']
    n_input=params['dims']
    dims=n_input
    chrom_list=list(range(2,23))#[19,20,21]#list(range(2,23)) #params['chrom'].split(':') 
    #chrom_list=list(range(int(chrom_list[0]),int(chrom_list[1])+1))
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    weights,biases,tensors=get_tensors([5,128,2],learning_rate)
    (x0, x1,x2,gt, accuracy, cost, optimizer, prob, rate)=tensors

    if params['val']:
        val_list=[]
        for v_path in params['test_path'].split(':'):
            
            if input_type=='path':
                _, x0_train, x1_train, x2_train, train_gt= get_data(params['test_path']+'pos','train',cpu=cpu)
                _, nx0_train, nx1_train, nx2_train, ntrain_gt= get_data(params['test_path']+'neg.high','train',cpu=cpu)
                
                vx0_test,vx1_test=None,None
                vx0_test=np.vstack([x0_train, nx0_train])
                vx1_test=np.vstack([x1_train, nx1_train])
                vx2_test=np.vstack([x2_train, nx2_train])
                vtest_gt=np.vstack([train_gt, ntrain_gt])
            
            else:
                vx0_test, vx1_test,vx2_test, vtest_gt=data['test']
            
            vx0_test=vx0_test.astype(np.float32)
            vx0_test[:,:,:,0]=vx0_test[:,:,:,0]-np.sum(vx0_test[:,:,:,0],axis=1)[:,np.newaxis,:]*vx0_test[:,:,:,1]
            
            vx1_test=vx1_test.astype(np.float32)
            vx1_test[:,:,:,0]=vx1_test[:,:,:,0]-np.sum(vx1_test[:,:,:,0],axis=1)[:,np.newaxis,:]*vx1_test[:,:,:,1]
            
            vx2_test=vx2_test.astype(np.float32)#[:,:,:64,:]
            vx2_test[:,:,:,0]=vx2_test[:,:,:,0]-np.sum(vx2_test[:,:,:,0],axis=1)[:,np.newaxis,:]*vx2_test[:,:,:,1]
            #vx2_test=vx2_test[:,:,:,0][:,:,:,np.newaxis]
            val_list.append((v_path,(vx0_test, vx1_test,vx2_test, vtest_gt)))
        
    
    init = tf.compat.v1.global_variables_initializer()
    saver = tf.compat.v1.train.Saver(max_to_keep=100)
        
    
    n_size=1
    with tf.compat.v1.Session(config=config)  as sess:
        sess.run(init)
        sess.run( tf.compat.v1.local_variables_initializer())
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
                f_path=os.path.join(params['train_path'],'chr%d/chr%d.pileups_nonrep.' %(chrom,chrom))
                
                #pos,mat,gt,allele_1,allele_2,ref
                    
                if input_type=='path':
                    _, x0_train, x1_train, x2_train, train_gt= get_data(f_path+'pos','train',cpu=cpu)

                    _, nx0_train, nx1_train, nx2_train, ntrain_gt=get_data(f_path+'neg.combined','train',cpu=cpu)
                
                else:
                    x0_train, x1_train, x2_train, train_gt=data['chr%d' %chrom]['pos']
                    
                    nx0_train, nx1_train, nx2_train, ntrain_gt=data['chr%d' %chrom]['neg']
                

                x0_train=x0_train.astype(np.float32)
                x0_train[:,:,:,0]=x0_train[:,:,:,0]-np.sum(x0_train[:,:,:,0],axis=1)[:,np.newaxis,:]*x0_train[:,:,:,1]
                
                x1_train=x1_train.astype(np.float32)
                x1_train[:,:,:,0]=x1_train[:,:,:,0]-np.sum(x1_train[:,:,:,0],axis=1)[:,np.newaxis,:]*x1_train[:,:,:,1]
                
                x2_train=x2_train.astype(np.float32)#[:,:,:64,:]
                x2_train[:,:,:,0]=x2_train[:,:,:,0]-np.sum(x2_train[:,:,:,0],axis=1)[:,np.newaxis,:]*x2_train[:,:,:,1]
                #x2_train=x2_train[:,:,:,0][:,:,:,np.newaxis]
                
                perm=np.random.permutation(len(nx2_train))
                
                nx0_train=np.take(nx0_train,perm,axis=0)
                nx1_train=np.take(nx1_train,perm,axis=0)
                
                nx2_train=np.take(nx2_train,perm,axis=0)
                ntrain_gt=np.take(ntrain_gt,perm,axis=0)
                
                nx0_train=nx0_train.astype(np.float32)
                nx0_train[:,:,:,0]=nx0_train[:,:,:,0]-np.sum(nx0_train[:,:,:,0],axis=1)[:,np.newaxis,:]*nx0_train[:,:,:,1]
                
                nx1_train=nx1_train.astype(np.float32)
                nx1_train[:,:,:,0]=nx1_train[:,:,:,0]-np.sum(nx1_train[:,:,:,0],axis=1)[:,np.newaxis,:]*nx1_train[:,:,:,1]
                
                nx2_train=nx2_train.astype(np.float32)#[:,:,:64,:]
                nx2_train[:,:,:,0]=nx2_train[:,:,:,0]-np.sum(nx2_train[:,:,:,0],axis=1)[:,np.newaxis,:]*nx2_train[:,:,:,1]
                #nx2_train=nx2_train[:,:,:,0][:,:,:,np.newaxis]
                
                chrom_data['chr%d'%chrom]=((x0_train,x1_train,x2_train, train_gt),(nx0_train, nx1_train, nx2_train, ntrain_gt))
                print('chromosome %d done | '%chrom,end='',flush=True)

        print('Finished reading pileups',flush=True)
        
        for k in range(iter_steps):
            
            for chrom in chrom_list:
                print('Training on chrom %d ' %(chrom),end='',flush=True)
                (x0_train, x1_train,x2_train, train_gt), (nx0_train, nx1_train, nx2_train, ntrain_gt) =chrom_data['chr%d'%chrom]


                n_start=-false_mltplr*len(x2_train)
                n_end=0
                
                training_loss=[]
                
                for i in range(iters):


                    n_start+=false_mltplr*len(x2_train)
                    n_end=n_start+false_mltplr*len(x2_train)

                    if n_end>len(nx2_train):
                        batch_nx0_train= np.vstack([nx0_train[n_start:,:,:,:],nx0_train[:n_end-len(nx0_train),:,:,:]])
                        batch_nx1_train= np.vstack([nx1_train[n_start:,:,:,:],nx1_train[:n_end-len(nx0_train),:,:,:]])
                        batch_nx2_train= np.vstack([nx2_train[n_start:,:,:,:],nx2_train[:n_end-len(nx2_train),:,:,:]])

                        
                        batch_ntrain_gt=np.vstack([ntrain_gt[n_start:,:],ntrain_gt[:n_end-len(nx2_train),:]])
                        

                        n_start=n_end-len(nx2_train)
                        n_end=n_start+false_mltplr*len(x2_train)
                    else:    
                        batch_nx0_train,batch_nx1_train,batch_nx2_train,batch_ntrain_gt= \
                        nx0_train[n_start:n_end,:,:,:],nx1_train[n_start:n_end,:,:,:],nx2_train[n_start:n_end,:,:,:],\
                        ntrain_gt[n_start:n_end,:]
                        
                        batch_nx2_train,batch_ntrain_gt= nx2_train[n_start:n_end,:,:,:], ntrain_gt[n_start:n_end,:]


                    total_false_size=int(false_mltplr*batch_size)

                    for batch in range(int(np.ceil(len(x2_train)/batch_size))):
                        batch_x0 = np.vstack([x0_train[batch*batch_size:min((batch+1)*batch_size,len(x0_train))],\
                                  batch_nx0_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx0_train))]])

                        batch_x1 = np.vstack([x1_train[batch*batch_size:min((batch+1)*batch_size,len(x1_train))],\
                                  batch_nx1_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx1_train))]])
                        
                        batch_x2 = np.vstack([x2_train[batch*batch_size:min((batch+1)*batch_size,len(x2_train))],\
                                  batch_nx2_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx2_train))]])
                        
                        batch_x=np.hstack([batch_x0, batch_x1, batch_x2])
                        
                        batch_gt = np.vstack([train_gt[batch*batch_size :min((batch+1)*batch_size,\
                                                   len(train_gt))], batch_ntrain_gt[total_false_size*batch : \
                                                   min(total_false_size*(batch+1), len(batch_ntrain_gt))]])
                        
                        
                       
                        opt,loss = sess.run([optimizer,cost], feed_dict={x2: batch_x, gt:batch_gt, rate:0.5})
                        
                        #training_loss.append(loss*len(batch_x))
                    
                    print('.',end='',flush=True)
                
                
                if (k<3 or chrom==chrom_list[-1]):
                    train_loss, train_acc, total_train_num=0, 0, 0
                    for batch in range(int(np.ceil(len(x2_train)/batch_size))):
                            batch_x0 = np.vstack([x0_train[batch*batch_size:min((batch+1)*batch_size,len(x0_train))],\
                                      batch_nx0_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                      len(batch_nx0_train))]])

                            batch_x1 = np.vstack([x1_train[batch*batch_size:min((batch+1)*batch_size,len(x1_train))],\
                                      batch_nx1_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                      len(batch_nx1_train))]])

                            batch_x2 = np.vstack([x2_train[batch*batch_size:min((batch+1)*batch_size,len(x2_train))],\
                                  batch_nx2_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx2_train))]])
                            
                            batch_x=np.hstack([batch_x0, batch_x1, batch_x2])
                        
                            batch_gt = np.vstack([train_gt[batch*batch_size :min((batch+1)*batch_size,\
                                                   len(train_gt))], batch_ntrain_gt[total_false_size*batch : \
                                                   min(total_false_size*(batch+1), len(batch_ntrain_gt))]])

                            batch_loss,batch_acc= sess.run([cost, accuracy], feed_dict={x2: batch_x, gt:batch_gt, rate:0.0})
                            total_train_num+=len(batch_x2)
                            train_loss+=batch_loss*len(batch_x2)
                            train_acc+=batch_acc

                    print('\n\n Training Loss= %.4f\n Training Accuracy= %.4f' %(train_loss/total_train_num, train_acc/total_train_num))
                
                
                if params['val'] and (k<3 or chrom==chrom_list[-1]):

                    for val in val_list:
                        print('\n')
                        print(val[0])
                        
                        vx0_test, vx1_test,vx2_test, vtest_gt=val[1]


                        v_loss,v_acc,total_vtest_num=0, 0, 0


                        for batch in range(int(np.ceil(len(vx2_test)/(batch_size)))):
                            vbatch_x0 = vx0_test[batch*batch_size:min((batch+1)*batch_size,len(vx0_test))]
                            vbatch_x1 = vx1_test[batch*batch_size:min((batch+1)*batch_size,len(vx1_test))]
                            vbatch_x2 = vx2_test[batch*batch_size:min((batch+1)*batch_size,len(vx2_test))]
                            
                            vbatch_x=np.hstack([vbatch_x0, vbatch_x1, vbatch_x2])
                            
                            vbatch_gt = vtest_gt[batch*batch_size:min((batch+1)*batch_size,len(vx2_test))]
                            
                            batch_loss,batch_acc= sess.run([cost, accuracy], feed_dict={x2: vbatch_x, gt: vbatch_gt, rate:0.0})

                            v_loss+=batch_loss*len(vbatch_x2)
                            v_acc+=batch_acc
                            total_vtest_num+=len(vbatch_x2)
                            
                        print('valid loss= %.4f\nvalid accuracy= %.4f \n' %(v_loss/total_vtest_num, v_acc/total_vtest_num))
                       
                        print(60*'-'+'\n')


            saver.save(sess, save_path=params['model'],global_step=save_num)
            elapsed=time.time()-t
            
            print ('Time Taken for Iteration %d-%d: %.2f seconds\n'\
                   %((save_num-1)*iters,save_num*iters,elapsed), flush=True)
            
            save_num+=1
            t=time.time()
            
        
        

def test_model(params,suffix='',prob_save=False):
    
    rev_allele_map={0:'N',1:'D',2:'I'}
    
    model_path,test_path,n_input,chrom,vcf_path= params['model'], params['test_path'],params['dims'],params['chrom'],params['vcf_path']
    
    tf.reset_default_graph()
    
    dims=n_input[:]

    weights,biases,tensors=get_tensors([5,128,2],0.0)
    (x0, x1,x2,gt, accuracy, cost, optimizer, prob, rate)=tensors
    
    init = tf.compat.v1.global_variables_initializer()
    sess = tf.compat.v1.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    

    batch_size=100
    total=[]
    rec_size=15371
    chnk=4
    neg_file=open(vcf_path+'.stats','w')
    neg_file.write('pos,ref\n')
    vcf_path=vcf_path
    
    with open(vcf_path,'w') as f:
        f_path=params['test_path']
        
        statinfo = os.stat(f_path)
        sz=statinfo.st_size        
        
        tmp_sz=list(range(0,sz,rec_size*32000))
        tmp_sz=tmp_sz[:chnk]
        tmp_sz=tmp_sz+[sz] if tmp_sz[-1]!=sz else tmp_sz
        total_prob=[]
        total_gt_prob=[]
        total_ref_=[]
        
        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        
        
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Probability">\n')
        f.write('##FORMAT=<ID=TP,Number=1,Type=String,Description="Indel type">\n')
        f.write('##gitcommit:%s\n' %str(git.Repo("/home/ahsanm/repos/NanoVar").heads[0].commit))
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
        
        prev=0
        prev_len=0
        
        for i in range(len(tmp_sz)-1):
            
            pos, x0_test, x1_test, x2_test= get_data(f_path,'test',a=tmp_sz[i], b=tmp_sz[i+1],cpu=params['cpu'])

            for batch in range(int(np.ceil(len(x0_test)/batch_size))):
                batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]

                batch_x0 = x0_test[batch*batch_size:min((batch+1)*batch_size,len(x0_test))]
                batch_x1 = x1_test[batch*batch_size:min((batch+1)*batch_size,len(x1_test))]
                batch_x2 = x2_test[batch*batch_size:min((batch+1)*batch_size,len(x2_test))]

                batch_x0=batch_x0.astype(np.float32)
                batch_x0[:,:,:,0]=batch_x0[:,:,:,0]/(np.sum(batch_x0[:,:,:,0],axis=1)[:,np.newaxis,:])-batch_x0[:,:,:,1]

                batch_x1=batch_x1.astype(np.float32)
                batch_x1[:,:,:,0]=batch_x1[:,:,:,0]/(np.sum(batch_x1[:,:,:,0],axis=1)[:,np.newaxis,:])-batch_x1[:,:,:,1]

                batch_x2=batch_x2.astype(np.float32)
                batch_x2[:,:,:,0]=batch_x2[:,:,:,0]/(np.sum(batch_x2[:,:,:,0],axis=1)[:,np.newaxis,:])-batch_x2[:,:,:,1]
                
                batch_x=np.hstack([batch_x0, batch_x1, batch_x2])
                batch_prob_all= sess.run(prob, feed_dict={x2: batch_x, rate:0.0})
                batch_pred_all=np.argmax(batch_prob_all,axis=1)
                
                
                '''batch_x=np.hstack([batch_x0, batch_x0, batch_x0])
                batch_prob_x0= sess.run(prob, feed_dict={x2: batch_x, rate:0.0})
                batch_pred_x0=np.argmax(batch_prob_x0,axis=1)
                
                batch_x=np.hstack([batch_x0, batch_x0, batch_x0])
                batch_prob_x1= sess.run(prob, feed_dict={x2: batch_x, rate:0.0})
                batch_pred_x1=np.argmax(batch_prob_x1,axis=1)'''
                
                qual_all=-10*np.log10(1e-6+1-batch_prob_all)
                
                for j in range(len(batch_pred_all)):
                    neg_file.write('%s,%d,%s,%.4f,%.4f,%.4f,%.4f\n' %(chrom,batch_pos[j], rev_gt_map[batch_pred_all[j]], batch_prob_all[j,0], batch_prob_all[j,1], batch_prob_all[j,2], batch_prob_all[j,3]))
                    
                    q=min(999,qual_all[j,batch_pred_all[j]])
                    if batch_pos[j]>prev:
                        
                        if batch_pred_all[j]==0 and batch_prob_all[j,0]<=0.95:
                            q=-10*np.log10(1e-6+batch_prob_all[j,0])
                            allele0_data=allele_prediction(batch_pos[j],batch_x0[j])
                            allele1_data=allele_prediction(batch_pos[j],batch_x1[j])

                            if allele0_data[0] and allele1_data[0]:
                                if allele0_data[0]==allele1_data[0] and allele0_data[1]==allele1_data[1]:
                                    s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP\t1/1:%s\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1],q,rev_gt_map[batch_pred_all[j]])
                                    f.write(s)
                                    prev=batch_pos[j]+max(len(allele0_data[0]), len(allele0_data[1]))

                                else:
                                    ref1,alt1=allele0_data
                                    ref2,alt2=allele1_data

                                    l=min(len(ref1),len(ref2))

                                    if len(ref1)>len(ref2):
                                        ref=ref1
                                        alt2=alt2+ref1[l:]

                                    else:
                                        ref=ref2
                                        alt1=alt1+ref2[l:]
                                    s='%s\t%d\t.\t%s\t%s,%s\t%.2f\tPASS\t.\tGT:TP\t1|2:%s\n' %(chrom, batch_pos[j], ref, alt1, alt2, q, rev_gt_map[batch_pred_all[j]])
                                    f.write(s)
                                    prev=batch_pos[j]+max(len(ref), len(alt1),len(alt2))

                            elif allele0_data[0]:
                                s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP\t0|1:%s\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], q, rev_gt_map[batch_pred_all[j]])
                                f.write(s)
                                prev=batch_pos[j]+max(len(allele0_data[0]), len(allele0_data[1]))

                            elif allele1_data[0]:
                                s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP\t1|0:%s\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1], q, rev_gt_map[batch_pred_all[j]])
                                f.write(s)
                                prev=batch_pos[j]+max(len(allele1_data[0]), len(allele1_data[1]))
                        
                        
                        elif batch_pred_all[j]>0:

                            if batch_pred_all[j]==1:

                                    allele_data=allele_prediction(batch_pos[j],batch_x2[j])
                                    if allele_data[0]:
                                        s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP\t1/1:%s\n' %(chrom, batch_pos[j], allele_data[0], allele_data[1], q, rev_gt_map[batch_pred_all[j]])
                                        f.write(s)
                                        prev=batch_pos[j]+max(len(allele_data[0]), len(allele_data[1]))

                            else: #(prediction 2 or 3)
                                allele0_data=allele_prediction(batch_pos[j],batch_x0[j])
                                allele1_data=allele_prediction(batch_pos[j],batch_x1[j])

                                if allele0_data[0] and allele1_data[0]: #(if two alleles predicted)
                                    
                                    if allele0_data[0]==allele1_data[0] and allele0_data[1]==allele1_data[1]: #(two alleles are same)
                                        s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP\t1/1:%s\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1],q,rev_gt_map[batch_pred_all[j]])
                                        f.write(s)
                                        prev=batch_pos[j]+max(len(allele0_data[0]), len(allele0_data[1]))

                                    else:
                                        ref1,alt1=allele0_data
                                        ref2,alt2=allele1_data

                                        l=min(len(ref1),len(ref2))

                                        if len(ref1)>len(ref2):
                                            ref=ref1
                                            alt2=alt2+ref1[l:]

                                        else:
                                            ref=ref2
                                            alt1=alt1+ref2[l:]
                                        s='%s\t%d\t.\t%s\t%s,%s\t%.2f\tPASS\t.\tGT:TP\t1|2:%s\n' %(chrom, batch_pos[j], ref, alt1, alt2, q, rev_gt_map[batch_pred_all[j]])
                                        f.write(s)
                                        prev=batch_pos[j]+max(len(ref), len(alt1),len(alt2))

                                elif allele0_data[0]:
                                    s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP\t0|1:%s\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], q, rev_gt_map[batch_pred_all[j]])
                                    f.write(s)
                                    prev=batch_pos[j]+max(len(allele0_data[0]), len(allele0_data[1]))

                                elif allele1_data[0]:
                                    s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP\t1|0:%s\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1], q, rev_gt_map[batch_pred_all[j]])
                                    f.write(s)
                                    prev=batch_pos[j]+max(len(allele1_data[0]), len(allele1_data[1]))
                batch_x0,batch_x1,batch_x2,batch_x=None,None,None,None
            pos, x0_test, x1_test, x2_test=None,None,None,None
    neg_file.close()
    
    
def allele_prediction(pos,mat):
    tmp_mat=mat[:,:,0]+mat[:,:,1]
    tmp_mat[4,:]=tmp_mat[4,:]-0.1
    ref_seq=''.join([rev_base_map[x] for x in np.argmax(mat[:,:,1],axis=0)])
    
    alt=''.join([rev_base_map[x] for x in np.argmax(tmp_mat,axis=0)])
    
    res=pairwise(alt.replace('-',''),ref_seq.replace('-',''))
    try:
        
        alt_new=res[0][0]
        ref_new=res[0][1]

        for i in range(60):
            if ref_new[i:i+5]==alt_new[i:i+5]:
                break

        i+=1
        s2=alt_new[:i].replace('-','')
        s1=ref_new[:i].replace('-','')
    
    
        s1=s1[0]+'.'+s1[1:]+'|'
        s2=s2[0]+'.'+s2[1:]+'|'
    except IndexError:
        return (None,None)
    if s1==s2:
        return (None,None)
    i=-1
    
    l=min(len(s1),len(s2))
    
    while s1[i]==s2[i] and s1[i]!='.' and s2[i]!='.':
        i-=1

    ref_out=s1[:i+1].replace('.','')
    allele_out=s2[:i+1].replace('.','')
    
    return (ref_out,allele_out) 

def read_pileups_from_file(options):
    rev_map={0:('I','I'), 1:('D','D'), 2:('N','I'), 3:('I','N'), 4:('N','D'), 5:('D','N'),6:('I','D'), 7:('D','I'), 8:('N','N')}
    allele_map={'N':0,'D':1,'I':2}
    rev_allele_map={0:'N',1:'D',2:'I'}
    
    fname,n,mode,dims=options
    file= open(fname,'r')
    file.seek(n)
    
    mat_0,mat_1,mat_2=[],[],[]
    pos=[]
    gt=[]
    
    i=0
    if mode=='train':

        while i<1000:
            i+=1
            c=file.read(12)
            if not c:
                break
            pos.append(int(c[:11]))

            
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

        mat_0=np.array(mat_0)
        mat_1=np.array(mat_1)
        mat_2=np.array(mat_2)
        pos=np.array(pos)
        
        gt=np.array(gt)
        #gt[gt==3]=2
        gt=np.eye(4)[gt].astype(np.int8)
        
        return (pos,mat_0,mat_1,mat_2,gt)

    else:
        
        while i<1000:
            i+=1
            c=file.read(11)
            if not c:
                break
            try:
                pos.append(int(c[:11]))

                m=file.read(dims[0]*dims[1]*dims[2]*4)
                p_mat_0=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
                mat_0.append(p_mat_0)

                m=file.read(dims[0]*dims[1]*dims[2]*4)
                p_mat_1=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
                mat_1.append(p_mat_1)

                m=file.read(dims[0]*dims[1]*dims[2]*4)
                p_mat_2=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
                mat_2.append(p_mat_2)
            except ValueError:
                print('error value for %s' %c )
        mat_0=np.array(mat_0)
        mat_1=np.array(mat_1)
        mat_2=np.array(mat_2)
        
        pos=np.array(pos)

        return (pos,mat_0,mat_1,mat_2)

    
def get_data(fname,mode,a=None, b=None,cpu=8):
    t=time.time()
    l=os.stat(fname).st_size
    dims=(5,128,2)
    
    mat_0,mat_1,gt=None,None,None
    
    if mode=='train':
        rec_size=12+dims[0]*dims[1]*dims[2]*4*3
        if a!=None and b!=None:
            my_array=[(fname,x,'train',dims) for x in range(a,b,1000*rec_size)]
        else:
            my_array=[(fname,x,'train',dims) for x in range(0,l,1000*rec_size)]
            
    else:
        rec_size=11+dims[0]*dims[1]*dims[2]*4*3
        if a!=None and b!=None:
            my_array=[(fname,x,'test',dims) for x in range(a,b,1000*rec_size)]
        else:
            my_array=[(fname,x,'test',dims) for x in range(0,l,1000*rec_size)]

    pool = mp.Pool(processes=cpu)
    results = pool.map(read_pileups_from_file, my_array)
    pool.close()  
    pool.join() 
    
    pos=np.vstack([res[0][:,np.newaxis] for res in results])
    mat_0=np.vstack([res[1] for res in results])
    mat_1=np.vstack([res[2] for res in results])
    mat_2=np.vstack([res[3] for res in results])
    
    if mode=='train':
        gt=np.vstack([res[4] for res in results])
        return pos,mat_0,mat_1,mat_2,gt
    else:
        return pos,mat_0,mat_1,mat_2
    
if __name__ == '__main__':
    print('git commit hash: %s' %str(git.Repo("/home/ahsanm/repos/NanoVar").heads[0].commit))
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
                'window':args.window,'ratio':args.ratio,'rt_path':args.rt_path,'cpu':args.cpu}
        genotype_caller_skinny(in_dict,neg_part=args.neg_part)
    
    else:
        in_dict={'cpu':args.cpu,'dims':input_dims,'test_path':args.test,'model':args.model, 'chrom':args.chrom, 'vcf_path':args.vcf, 'workdir':args.workdir}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
