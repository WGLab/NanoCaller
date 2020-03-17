import time,os,copy,argparse,subprocess,psutil
import pandas as pd
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect_indel_real import *
from matplotlib import pyplot as plt
import git

#config = tf.ConfigProto(device_count={"CPU": 64})

config = tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True

rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}
rev_map={0:('I','I'), 1:('D','D'), 2:('N','I'), 3:('I','N'), 4:('N','D'), 5:('D','N'),6:('I','D'), 7:('D','I'), 8:('N','N')}
allele_map={'N':0,'D':1,'I':2}
rev_allele_map={0:'N',1:'D',2:'I'}

def genotype_caller_skinny(params,input_type='path',data=None,attempt=0,neg_part='neg.combined'):
    tf.reset_default_graph()
    
    false_mltplr=1
    
    cpu=params['cpu']
    n_input=params['dims']
    dims=n_input
    chrom_list=list(range(2,23))#[19,20,21]#list(range(2,23)) #params['chrom'].split(':') 
    #chrom_list=list(range(int(chrom_list[0]),int(chrom_list[1])+1))
    
    training_iters, learning_rate, batch_size= params['iters'],\
    params['rate'], params['size']

    weights,biases,tensors=get_tensors(n_input,learning_rate)
    (x0, x1, label_1,label_2, accuracy, accuracy_0, accuracy_1, cost, optimizer, prob_0,prob_1, rate)=tensors

    if params['val']:
        val_list=[]
        for v_path in params['test_path'].split(':'):
            
            if input_type=='path':
                _, x0_train, x1_train, train_allele_1,train_allele_2= read_pileups_from_file((params['test_path']+'pos',0, 'train', dims))
                _, nx0_train, nx1_train, ntrain_allele_1,ntrain_allele_2=read_pileups_from_file((params['test_path']+'neg.high',0, 'train', dims))

                vx0_test=np.vstack([x0_train, nx0_train])
                vx1_test=np.vstack([x1_train, nx1_train])
                vtest_allele_1=np.vstack([train_allele_1, ntrain_allele_1])
                vtest_allele_2=np.vstack([train_allele_2, ntrain_allele_2])
            
            else:
                vx0_test, vx1_test, vtest_allele_1,vtest_allele_2=data['test']
            
            vx0_test=vx0_test.astype(np.float32)
            vx0_test[:,:,:,0]=vx0_test[:,:,:,0]/(np.sum(vx0_test[:,:,:,0],axis=1)[:,np.newaxis,:])-vx0_test[:,:,:,1]
            
            vx1_test=vx1_test.astype(np.float32)
            vx1_test[:,:,:,0]=vx1_test[:,:,:,0]/(np.sum(vx1_test[:,:,:,0],axis=1)[:,np.newaxis,:])-vx1_test[:,:,:,1]
                
            val_list.append((v_path,(vx0_test, vx1_test, vtest_allele_1,vtest_allele_2)))
        
    
    init = tf.compat.v1.global_variables_initializer()
    saver = tf.compat.v1.train.Saver(max_to_keep=100)
        
    
    n_size=1
    with tf.compat.v1.Session(config=config)  as sess:
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
                f_path=os.path.join(params['train_path'],'chr%d/chr%d.pileups_mixed.' %(chrom,chrom))
                
                #pos,mat,gt,allele_1,allele_2,ref
                    
                if input_type=='path':
                    _, x0_train, x1_train, train_allele_1,train_allele_2= read_pileups_from_file((f_path+'pos',0,'train', n_input))

                    _, nx0_train, nx1_train, ntrain_allele_1,ntrain_allele_2=read_pileups_from_file((f_path+'neg.combined', 0,'train', n_input))
                
                else:
                    x0_train, x1_train, train_allele_1,train_allele_2=data['chr%d' %chrom]['pos']
                    
                    #x0_train=x0_train
                    #, x1_train
                    
                    
                    nx0_train, nx1_train, ntrain_allele_1,ntrain_allele_2=data['chr%d' %chrom]['neg']
                

                x0_train=x0_train.astype(np.float32)
                x0_train[:,:,:,0]=x0_train[:,:,:,0]/(np.sum(x0_train[:,:,:,0],axis=1)[:,np.newaxis,:])-x0_train[:,:,:,1]
                
                x1_train=x1_train.astype(np.float32)
                x1_train[:,:,:,0]=x1_train[:,:,:,0]/(np.sum(x1_train[:,:,:,0],axis=1)[:,np.newaxis,:])-x1_train[:,:,:,1]
                
                perm=np.random.permutation(len(nx0_train))
                nx0_train=np.take(nx0_train,perm,axis=0)
                nx1_train=np.take(nx1_train,perm,axis=0)
                ntrain_allele_1=np.take(ntrain_allele_1,perm,axis=0)
                ntrain_allele_2=np.take(ntrain_allele_2,perm,axis=0)
                
                nx0_train=nx0_train.astype(np.float32)
                nx0_train[:,:,:,0]=nx0_train[:,:,:,0]/(np.sum(nx0_train[:,:,:,0],axis=1)[:,np.newaxis,:])-nx0_train[:,:,:,1]
                
                nx1_train=nx1_train.astype(np.float32)
                nx1_train[:,:,:,0]=nx1_train[:,:,:,0]/(np.sum(nx1_train[:,:,:,0],axis=1)[:,np.newaxis,:])-nx1_train[:,:,:,1]
                
                chrom_data['chr%d'%chrom]=((x0_train, x1_train, train_allele_1,train_allele_2),(nx0_train, nx1_train, ntrain_allele_1,ntrain_allele_2))
                print('chromosome %d done | '%chrom,end='',flush=True)

        print('Finished reading pileups',flush=True)
        
        for k in range(iter_steps):
            
            for chrom in chrom_list:
                print('Training on chrom %d ' %(chrom),end='',flush=True)
                (x0_train, x1_train, train_allele_1,train_allele_2), (nx0_train, nx1_train, ntrain_allele_1,ntrain_allele_2)=chrom_data['chr%d'%chrom]


                n_start=-false_mltplr*len(x0_train)
                n_end=0
                
                training_loss=[]
                
                for i in range(iters):


                    n_start+=false_mltplr*len(x0_train)
                    n_end=n_start+false_mltplr*len(x0_train)

                    if n_end>len(nx0_train):
                        batch_nx0_train= np.vstack([nx0_train[n_start:,:,:,:],nx0_train[:n_end-len(nx0_train),:,:,:]])
                        batch_nx1_train= np.vstack([nx1_train[n_start:,:,:,:],nx1_train[:n_end-len(nx0_train),:,:,:]])

                        batch_ntrain_allele_1=np.vstack([ntrain_allele_1[n_start:,:],ntrain_allele_1[:n_end-len(nx0_train),:]])
                        batch_ntrain_allele_2=np.vstack([ntrain_allele_2[n_start:,:],ntrain_allele_2[:n_end-len(nx0_train),:]])
                        

                        n_start=n_end-len(nx0_train)
                        n_end=n_start+false_mltplr*len(x0_train)
                    else:    
                        batch_nx0_train,batch_nx1_train,batch_ntrain_allele_1, batch_ntrain_allele_2= \
                        nx0_train[n_start:n_end,:,:,:],nx1_train[n_start:n_end,:,:,:],\
                        ntrain_allele_1[n_start:n_end,:],ntrain_allele_2[n_start:n_end,:]


                    total_false_size=int(false_mltplr*batch_size)

                    for batch in range(int(np.ceil(len(x0_train)/batch_size))):
                        batch_x0 = np.vstack([x0_train[batch*batch_size:min((batch+1)*batch_size,len(x0_train))],\
                                  batch_nx0_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx0_train))]])

                        batch_x1 = np.vstack([x1_train[batch*batch_size:min((batch+1)*batch_size,len(x1_train))],\
                                  batch_nx1_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                  len(batch_nx1_train))]])
                        
                        batch_allele_1 = np.vstack([train_allele_1[batch*batch_size :min((batch+1)*batch_size,\
                                                   len(train_allele_1))], batch_ntrain_allele_1[total_false_size*batch : \
                                                   min(total_false_size*(batch+1), len(batch_ntrain_allele_1))]])
                        
                        batch_allele_2 = np.vstack([train_allele_2[batch*batch_size :min((batch+1)*batch_size,\
                                                   len(train_allele_2))], batch_ntrain_allele_2[total_false_size*batch : \
                                                   min(total_false_size*(batch+1), len(batch_ntrain_allele_2))]])

                       
                        opt,loss = sess.run([optimizer,cost], feed_dict={x0: batch_x0, x1: batch_x1, label_1:batch_allele_1, label_2:batch_allele_2, rate:0.5})
                        
                        #training_loss.append(loss*len(batch_x))
                    
                    print('.',end='',flush=True)
                
                
                if (k<3 or chrom==chrom_list[-1]):
                    train_loss, train_acc, total_train_num=0, 0, 0
                    for batch in range(int(np.ceil(len(x0_train)/batch_size))):
                            batch_x0 = np.vstack([x0_train[batch*batch_size:min((batch+1)*batch_size,len(x0_train))],\
                                      batch_nx0_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                      len(batch_nx0_train))]])

                            batch_x1 = np.vstack([x1_train[batch*batch_size:min((batch+1)*batch_size,len(x1_train))],\
                                      batch_nx1_train[ total_false_size*batch : min(total_false_size*(batch+1),\
                                      len(batch_nx1_train))]])

                            batch_allele_1 = np.vstack([train_allele_1[batch*batch_size :min((batch+1)*batch_size,\
                                                       len(train_allele_1))], batch_ntrain_allele_1[total_false_size*batch : \
                                                       min(total_false_size*(batch+1), len(batch_ntrain_allele_1))]])

                            batch_allele_2 = np.vstack([train_allele_2[batch*batch_size :min((batch+1)*batch_size,\
                                                       len(train_allele_2))], batch_ntrain_allele_2[total_false_size*batch : \
                                                       min(total_false_size*(batch+1), len(batch_ntrain_allele_2))]])

                            batch_loss,batch_acc= sess.run([cost, accuracy], feed_dict={x0: batch_x0, x1: batch_x1, label_1:batch_allele_1, label_2:batch_allele_2, rate:0.0})
                            total_train_num+=len(batch_x0)
                            train_loss+=batch_loss*len(batch_x0)
                            train_acc+=batch_acc

                    print('\n\n Training Loss= %.4f\n Training Accuracy= %.4f' %(train_loss/total_train_num, train_acc/total_train_num))
                
                
                if params['val'] and (k<3 or chrom==chrom_list[-1]):

                    for val in val_list:
                        print('\n')
                        print(val[0])
                        
                        vx0_test, vx1_test, vtest_allele_1,vtest_allele_2=val[1]


                        v_loss,v_acc,total_vtest_num=0, 0, 0


                        for batch in range(int(np.ceil(len(vx0_test)/(batch_size)))):
                            vbatch_x0 = vx0_test[batch*batch_size:min((batch+1)*batch_size,len(vx0_test))]
                            vbatch_x1 = vx1_test[batch*batch_size:min((batch+1)*batch_size,len(vx1_test))]
                            vbatch_allele_1 = vtest_allele_1[batch*batch_size:min((batch+1)*batch_size,len(vx0_test))]
                            vbatch_allele_2 = vtest_allele_2[batch*batch_size:min((batch+1)*batch_size,len(vx0_test))]

                            batch_loss,batch_acc= sess.run([cost, accuracy], feed_dict={x0: vbatch_x0, x1: vbatch_x1, label_1: vbatch_allele_1, label_2: vbatch_allele_2, rate:0.0})

                            v_loss+=batch_loss*len(vbatch_x0)
                            v_acc+=batch_acc
                            total_vtest_num+=len(vbatch_x0)
                            
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

    weights,biases,tensors=get_tensors(n_input,0.0)
    (x0, x1, label_1,label_2, accuracy, accuracy_0, accuracy_1, cost, optimizer, prob_0,prob_1, rate)=tensors
    
    init = tf.compat.v1.global_variables_initializer()
    sess = tf.compat.v1.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    

    batch_size=100
    total=[]
    rec_size=5131
    chnk=1
    neg_file=open(vcf_path+'.stats','w')
    neg_file.write('pos,ref\n')
    vcf_path=vcf_path
    
    with open(vcf_path,'w') as f:
        f_path=params['test_path']
        
        statinfo = os.stat(f_path)
        sz=statinfo.st_size        
        
        tmp_sz=list(range(0,sz,rec_size*(sz//(chnk*rec_size))))
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
        f.write('##gitcommit:%s\n' %str(git.Repo("/home/ahsanm/repos/NanoVar").heads[0].commit))
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
        
        
        for i in range(len(tmp_sz)-1):
            
            pos, x0_test, x1_test= get_data(f_path,a=tmp_sz[i], b=tmp_sz[i+1])
            

            x0_test=x0_test.astype(np.float32)
            x0_test[:,:,:,0]=x0_test[:,:,:,0]/(np.sum(x0_test[:,:,:,0],axis=1)[:,np.newaxis,:])-x0_test[:,:,:,1]

            x1_test=x1_test.astype(np.float32)
            x1_test[:,:,:,0]=x1_test[:,:,:,0]/(np.sum(x1_test[:,:,:,0],axis=1)[:,np.newaxis,:])-x1_test[:,:,:,1]

            for batch in range(int(np.ceil(len(x0_test)/batch_size))):
                batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]

                batch_x0 = x0_test[batch*batch_size:min((batch+1)*batch_size,len(x0_test))]
                batch_x1 = x1_test[batch*batch_size:min((batch+1)*batch_size,len(x1_test))]

                batch_prob_0,batch_prob_1= sess.run([prob_0,prob_1], feed_dict={x0: batch_x0, x1: batch_x1, rate:0.0})

                batch_pred_0=np.argmax(batch_prob_0,axis=1)
                batch_pred_1=np.argmax(batch_prob_1,axis=1)

                qual_0=-10*np.log10(1e-10+1-batch_prob_0)
                qual_1=-10*np.log10(1e-10+1-batch_prob_1)

                for j in range(len(batch_pred_0)):

                    if batch_pred_0[j]!=0 or batch_pred_1[j]!=0:
                        q0=min(99,qual_0[j,batch_pred_0[j]])
                        q1=min(99,qual_1[j,batch_pred_1[j]])
                        
                        allele0_data=allele_prediction(batch_pos[j],batch_x0[j],rev_allele_map[batch_pred_0[j]],30)
                            
                        allele1_data=allele_prediction(batch_pos[j],batch_x1[j],rev_allele_map[batch_pred_1[j]],30)
                        
                        if batch_pred_0[j]!=0 and allele0_data[0] and batch_pred_1[j]!=0 and allele1_data[0]:
                            q3=(q0+q1)/2
                            #s='%s,%d,%s,%s,%.4f,%.4f,%.4f,%.4f\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1],q0,q1,(q0+q1)/2,q3)
                            if allele0_data[1]==allele1_data[1]:
                                s='%s\t%d\t.\t%s\t%s\t%.2f\t%s\t.\tGT\t%s\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1],q3,'PASS','1/1' )
                            else:
                                s='%s\t%d\t.\t%s\t%s,%s\t%.2f\t%s\t.\tGT\t%s\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], allele1_data[1],q3,'PASS','1|2' )
                                
                            f.write(s)
                            
                            
                            
                        elif batch_pred_0[j]!=0 and allele0_data[0]:
                            q3=q0
                            
                            #s='%s,%d,%s,%s,%.4f,%.4f,%.4f,%.4f\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1],q0,q1,(q0+q1)/2,q3)
                            s='%s\t%d\t.\t%s\t%s\t%.2f\t%s\t.\tGT\t%s\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1],q3,'PASS','1|0' )
                            f.write(s)
                            
                        elif batch_pred_1[j]!=0 and allele1_data[0]:
                            q3=q1
                            
                            #s='%s,%d,%s,%s,%.4f,%.4f,%.4f,%.4f\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1],q0,q1,(q0+q1)/2,q3)
                            s='%s\t%d\t.\t%s\t%s\t%.2f\t%s\t.\tGT\t%s\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1], q3,'PASS','0|1' )
                            f.write(s)
                       
                    neg_file.write('%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n' %(batch_pos[j], batch_prob_0[j,0], batch_prob_0[j,1], batch_prob_0[j,2], batch_prob_1[j,0],batch_prob_1[j,1],batch_prob_1[j,2]))

def allele_prediction(pos,mat1,var_type,l):
    tmp_mat=mat1[:,:,0]+mat1[:,:,1]
    tmp_mat[4,:]=tmp_mat[4,:]-mat1[4,:,1]
    ref_mat=np.argmax(mat1[:,:,1],axis=0)
    
    res=[]
    
    if var_type=='N':
        return (None,None)
    
    elif var_type=='D':
        allele=''
        ref_allele=''

        for i in range(64):
            if ref_mat[i]!=4:
                ref_allele+=rev_base_map[ref_mat[i]]
                
            if max(abs(tmp_mat[:,i]))<0.5 or ref_mat[i]==4:
                continue
            else:                
                if max(tmp_mat[:4,i])>0.5:
                        allele+=rev_base_map[np.argmax(tmp_mat[:4,i])]
                        
            
            res.append((ref_allele,allele))
            
            if max(len(ref_allele),len(allele))>l:
                break
            
                    
    else:
        allele=''
        ref_allele=''
                        
        for i in range(64):
            if ref_mat[i]!=4:
                ref_allele+=rev_base_map[ref_mat[i]]
                
            if max(abs(tmp_mat[:,i]))<0.5:
                continue
            else:
                if ref_mat[i]==4 and max(tmp_mat[:4,i])>=0.5:
                    allele+=rev_base_map[np.argmax(tmp_mat[:4,i])]
                
                elif ref_mat[i]!=4: 
                    
                    if max(tmp_mat[:4,i])>=0.5:
                        allele+=rev_base_map[np.argmax(tmp_mat[:4,i])]
                    else:
                        allele+=rev_base_map[ref_mat[i]]
                        
            res.append((ref_allele,allele))  
            
            if max(len(ref_allele),len(allele))>l:
                break
    
    output=[]
    
    for pair in res:
        s1,s2=pair

        i=-1

        try:
            while s1[i]==s2[i]:
                i-=1
            if i==-1:
                ref_out=s1
                allele_out=s2
            else:
                ref_out=s1[:i+1]
                allele_out=s2[:i+1]
            output.append((ref_out,allele_out))
            
        except IndexError:
            pass
    if len(output)==0:
        return (None,None)
    
    else:
        if var_type=='I':
            output.sort(key=lambda x: len(x[0]))
            lst=[pair for pair in output if len(pair[0])==len(output[0][0])]
            lst.sort(key=lambda x: len(x[1]), reverse=True)

        elif var_type=='D':
            output.sort(key=lambda x: len(x[1]))
            lst=[pair for pair in output if len(pair[1])==len(output[0][1])]
            lst.sort(key=lambda x: len(x[0]), reverse=True)
            
        return lst[0]
    
    
    
def read_pileups_from_file(options):
    rev_map={0:('I','I'), 1:('D','D'), 2:('N','I'), 3:('I','N'), 4:('N','D'), 5:('D','N'),6:('I','D'), 7:('D','I'), 8:('N','N')}
    allele_map={'N':0,'D':1,'I':2}
    rev_allele_map={0:'N',1:'D',2:'I'}
    
    fname,n,mode,dims=options
    file= open(fname,'r')
    file.seek(n)
    
    mat_0,mat_1=[],[]
    pos=[]
    allele1,allele2=[],[]
    
    i=0
    if mode=='train':

        while True:
            i+=1
            c=file.read(12)
            if not c:
                break
            pos.append(int(c[:11]))

            
            allele1.append(allele_map[rev_map[int(c[11])][0]])
            allele2.append(allele_map[rev_map[int(c[11])][1]])
 
            m=file.read(dims[0]*dims[1]*dims[2]*4)
            p_mat_0=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
            mat_0.append(p_mat_0)
            
            m=file.read(dims[0]*dims[1]*dims[2]*4)
            p_mat_1=np.array([int(float(m[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1],dims[2]))
            mat_1.append(p_mat_1)

        mat_0=np.array(mat_0)
        mat_1=np.array(mat_1)
        pos=np.array(pos)
        allele1=np.eye(3)[np.array(allele1)].astype(np.int8)
        allele2=np.eye(3)[np.array(allele2)].astype(np.int8)
        
        return (pos,mat_0,mat_1,allele1,allele2)

    else:
        
        while i<1000:
            i+=1
            c=file.read(11)
            if not c:
                break
            
            try:
                pos.append(int(c[:11]))

                m0=file.read(dims[0]*dims[1]*dims[2]*4)
                p_mat_0=np.array([int(float(m0[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1], dims[2])).astype(np.int16)
                mat_0.append(p_mat_0)

                m1=file.read(dims[0]*dims[1]*dims[2]*4)
                p_mat_1=np.array([int(float(m1[4*i:4*i+4])) for i in range(dims[0]*dims[1]*dims[2])]).reshape((dims[0],dims[1], dims[2])).astype(np.int16)
                mat_1.append(p_mat_1)
            except ValueError:
                print(pos[-1],c,i)
                print(m0)
                print(m1)
                return
                
        mat_0=np.array(mat_0)
        mat_1=np.array(mat_1)
        pos=np.array(pos)

        return (pos,mat_0,mat_1)

    
def get_data(fname,a=None, b=None):
    t=time.time()
    l=os.stat(fname).st_size
    dims=(5,64,2)
    rec_size=5131
    if a!=None and b!=None:
        my_array=[(fname,x,'test',dims) for x in range(a,b,1000*rec_size)]
    else:
        my_array=[(fname,x,'test',dims) for x in range(0,l,1000*rec_size)]

    cpu=8
    pool = mp.Pool(processes=cpu)
    results = pool.map(read_pileups_from_file, my_array)
    pool.close()  
    pool.join() 
    
    pos=np.vstack([res[0][:,np.newaxis] for res in results])
    mat_0=np.vstack([res[1] for res in results])
    mat_1=np.vstack([res[2] for res in results])

            

    return pos,mat_0,mat_1
    
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
                'window':args.window,'ratio':args.ratio,'rt_path':args.rt_path}
        genotype_caller_skinny(in_dict,neg_part=args.neg_part)
    
    else:
        in_dict={'cpu':args.cpu,'dims':input_dims,'test_path':args.test,'model':args.model, 'chrom':args.chrom, 'vcf_path':args.vcf, 'workdir':args.workdir}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
