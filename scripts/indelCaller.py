import time,os,copy,argparse,subprocess
import pandas as pd
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect_indel import *
from Bio import pairwise2
from generate_indel_pileups import generate

config =  tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

rev_gt_map={0:'hom-ref', 1:'hom-alt', 2:'het-ref', 3:'het-alt'}
rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

def pairwise(x,y):
    alignments = pairwise2.align.localms(x, y, 2, -1.0, -0.5, -0.1)

    return alignments

def test_model(params,pool):
    
    rev_allele_map={0:'N',1:'D',2:'I'}
    
    vcf_path,prefix=params['vcf_path'], params['prefix']
    chrom,start,end=params['chrom'], params['start'],params['end']

    tf.reset_default_graph()

    weights,biases,tensors=get_tensors([5,128,2],0.0)
    (x0, x1,x2,gt, accuracy, cost, optimizer, prob, rate)=tensors
    
    init = tf.compat.v1.global_variables_initializer()
    sess = tf.compat.v1.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    dirname = os.path.dirname(__file__)
    saver.restore(sess, os.path.join(dirname, 'release_data/NanoCaller_indels/model-30'))
    
    batch_size=100

    neg_file=open(os.path.join(vcf_path,'%s.indel_stats' %prefix),'w')
    neg_file.write('pos,ref\n')
    
    with open(os.path.join(vcf_path,'%s.indels.vcf' %prefix),'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Probability">\n')
        f.write('##FORMAT=<ID=TP,Number=1,Type=String,Description="Indel type">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s\n' %params['sample'])
        
        prev=0
        prev_len=0
        
        for mbase in range(start,end,int(1e7)):
            d = copy.deepcopy(params)
            d['start']=mbase
            d['end']=min(end,mbase+int(1e7))
            
            pos, x0_test, x1_test, x2_test=generate(d,pool)
            
            if len(pos)==0:
                continue

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
            
            f.flush()
            os.fsync(f.fileno())
            neg_file.flush()
            os.fsync(neg_file.fileno())
            
    neg_file.close()
    outfile=os.path.join(vcf_path,'%s.indels' %prefix)
    stream=os.popen("bcftools sort %s.vcf|bgziptabix %s.vcf.gz" %(outfile, outfile))
    stream.read()
    return outfile
    
    
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
