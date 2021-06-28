import time,os,copy,argparse,subprocess, sys, pysam, datetime

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect import *
from generate_SNP_pileups import get_snp_testing_candidates
from intervaltree import Interval, IntervalTree
from utils import *

if type(tf.contrib) != type(tf): tf.contrib._warning = None
config = tf.ConfigProto()
config.gpu_options.allow_growth = True

num_to_base_map={0:'A',1:'G',2:'T',3:'C'}
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

snp_model_dict={'NanoCaller1':'release_data/ONT_models/SNPs/NanoCaller1_beta/model-rt-1',
                'NanoCaller2':'release_data/ONT_models/SNPs/NanoCaller1_beta/model-rt-1',
                'NanoCaller3':'release_data/clr_models/SNPs/NanoCaller3_beta/model-rt-100',
                'ONT-HG001':'release_data/ONT_models/SNPs/HG001_guppy4.2.2_giab-3.3.2/model-1',
                'ONT-HG001_GP2.3.8':'release_data/ONT_models/SNPs/HG001_guppy2.3.8_giab-3.3.2/model-100',
                'ONT-HG001_GP2.3.8-4.2.2':'release_data/ONT_models/SNPs/HG001_guppy2.3.8_guppy4.2.2_giab-3.3.2/model-100',
                'ONT-HG001-4_GP4.2.2':'release_data/ONT_models/SNPs/HG001_guppy4.2.2_giab-3.3.2_HG002-4_guppy4.2.2_giab-4.2.1/model-100',
                'ONT-HG002':'release_data/ONT_models/SNPs/HG002_guppy4.2.2_giab-4.2.1/model-100',
                'ONT-HG002_GP4.2.2_v3.3.2':'release_data/ONT_models/SNPs/HG002_guppy4.2.2_giab-3.3.2/model-100',
                'ONT-HG002_GP2.3.4_v3.3.2':'release_data/ONT_models/SNPs/HG002_guppy2.3.4_giab-3.3.2/model-100',
                'ONT-HG002_GP2.3.4_v4.2.1':'release_data/ONT_models/SNPs/HG002_guppy2.3.4_giab-4.2.1/model-100',
                'ONT-HG002_r10.3':'release_data/ONT_models/SNPs/HG002_r10.3_guppy4.0.11_giab-4.2.1/model-100',
                'ONT-HG002_bonito':'release_data/ONT_models/SNPs/HG002_bonito_giab-4.2.1/model-100',
                'CCS-HG001':'release_data/hifi_models/SNPs/HG001_giab-3.3.2/model-100',
                'CCS-HG002':'release_data/hifi_models/SNPs/HG002_giab-4.2.1/model-100',
                'CCS-HG001-4':'release_data/hifi_models/SNPs/HG001_giab-3.3.2_HG002-4_giab-4.2.1/model-100',
                'CLR-HG002':'release_data/clr_models/SNPs/HG002_giab-4.2.1/model-100'
               }

def get_SNP_model(snp_model):
    if os.path.exists(snp_model):
        if os.path.isdir(snp_model):
            snp_model_path=glob.glob(os.path.join(snp_model,'*.meta'))[0].rstrip('.meta')
    
    elif snp_model in snp_model_dict:
        dirname = os.path.dirname(__file__)
        snp_model_path=os.path.join(dirname, snp_model_dict[snp_model])
    
    else:
        return None,None
    
    coverage_path='%s.coverage' %snp_model_path
    
    if os.path.exists(coverage_path):
        train_coverage=float(open(coverage_path,'r').readlines()[0].rstrip('\n'))
    else:
        train_coverage=0
        
    return snp_model_path, train_coverage
        
def test_model(params,pool):
    print('%s: SNP calling started.' %(str(datetime.datetime.now())), flush=True)
    
    vcf_path,prefix= params['vcf_path'], params['prefix']
    chrom,start,end=params['chrom'], params['start'],params['end']          
    
    coverage=get_coverage(params, pool)
    
    print('%s: Coverage=%.2fx.' %(str(datetime.datetime.now()), coverage), flush=True)
    
    if coverage==0:
        print('%s: No coverage found for the contig %s.' %(str(datetime.datetime.now()), chrom), flush=True)
        return
    
    
    n_input=[5,41,5]

    tf.reset_default_graph()
    
    weights,biases,tensors=get_tensors(n_input,0.0)
    (x,GT_label,A_label, G_label, T_label, C_label,GT_score, A_score, G_score, T_score, C_score, accuracy_GT, accuracy_A,  accuracy_G,  accuracy_T,  accuracy_C, prediction_accuracy_GT, prediction_accuracy_A,  prediction_accuracy_G,  prediction_accuracy_T,  prediction_accuracy_C, prediction_GT, prediction_A,  prediction_G,  prediction_T,  prediction_C, accuracy, cost, optimizer,  cost_GT, cost_A, cost_G, cost_T, cost_C,A_ref,G_ref,T_ref,C_ref,prob_GT,prob_A,prob_G,prob_T,prob_C,keep)=tensors
    
    model_path, train_coverage=get_SNP_model(params['snp_model'])
    
    if model_path==None:
        print('Invalid SNP model name or path', flush=True)
        return
    
    train_coverage=coverage if train_coverage==0 else train_coverage
        
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)    
    
    batch_size=1000
    neg_file=open(os.path.join(vcf_path,'%s.snp_stats' %prefix),'w')
    neg_file.write('pos,ref,prob_GT,prob_A,prob_G,prob_T,prob_C,DP,freq\n')
    
    with open(os.path.join(vcf_path,'%s.snps.vcf' %prefix),'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        f.write('##FORMAT=<ID=FQ,Number=1,Type=Float,Description="Alternative Allele Frequency">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s\n' %params['sample'])
        
        in_dict_list=[]
        
        for mbase in range(start,end,200000):
            d = copy.deepcopy(params)
            d['start']=mbase
            d['end']=min(end,mbase+200000)
            in_dict_list.append(d)
        
        if params['cpu']==1:
            result=map(get_snp_testing_candidates, in_dict_list)
        
        else:
            result=pool.imap_unordered(get_snp_testing_candidates, in_dict_list)

        total_regions=len(in_dict_list)
        completed=0
        
        for res in result:
            
            pos,test_ref,x_test,dp,freq=res
            completed+=1
            
            if len(pos)==0:
                continue
       
            x_test=x_test.astype(np.float32)
            
            x_test[:,1:,:,:4]=x_test[:,1:,:,:4]*(train_coverage/coverage)

            for batch in range(int(np.ceil(len(x_test)/batch_size))):
                batch_freq=freq[batch*batch_size:min((batch+1)*batch_size,len(freq))]
                batch_dp=dp[batch*batch_size:min((batch+1)*batch_size,len(dp))]
                batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]

                batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]

                batch_prob_GT,batch_prob_A,batch_prob_G,batch_prob_C,batch_prob_T= sess.run([prob_GT,prob_A,prob_G,prob_C,prob_T],\
                                       feed_dict={x: batch_x, A_ref:batch_ref[:,0][:,np.newaxis], G_ref:batch_ref[:,1][:,np.newaxis], T_ref:batch_ref[:,2][:,np.newaxis], C_ref:batch_ref[:,3][:,np.newaxis],keep:1.0})

                batch_pred_GT=np.argmax(batch_prob_GT,axis=1)

                batch_probs=np.hstack([batch_prob_A[:,1][:,np.newaxis], batch_prob_G[:,1][:,np.newaxis], batch_prob_T[:,1][:,np.newaxis], batch_prob_C[:,1][:,np.newaxis]])


                batch_pred=np.argsort(batch_probs,axis=1)
                
                batch_ref_vec=batch_ref
                
                batch_ref=np.argmax(batch_ref,1)
                
                batch_pred_GT=np.sum(batch_probs>=0.5,axis=1)
                
                sort_probs=np.sort(batch_probs,axis=1)
                
                for j in range(len(batch_pred_GT)):

                    if batch_pred_GT[j]>=2: # if het
                            pred1,pred2=batch_pred[j,-1],batch_pred[j,-2]
                            if pred1==batch_ref[j]:
                                        s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t.\tGT:DP:FQ\t%s:%d:%.4f\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred2], min(999,-100*np.log10(1e-10+ 1-batch_probs[j,pred2])),'PASS','0/1', batch_dp[j], batch_freq[j])
                                        f.write(s)

                            elif pred2==batch_ref[j] and batch_probs[j,pred2]>=0.5:
                                s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t.\tGT:DP:FQ\t%s:%d:%.4f\n' %(chrom,batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred1], min(999,-100*np.log10(1e-10+ 1-batch_probs[j,pred2])),'PASS','1/0', batch_dp[j], batch_freq[j])
                                f.write(s)

                            elif pred2!=batch_ref[j] and pred1!=batch_ref[j] and batch_probs[j,pred2]>=0.5:
                                s='%s\t%d\t.\t%s\t%s,%s\t%.3f\t%s\t.\tGT:DP:FQ\t%s:%d:%.4f\n' %\
                    (chrom,batch_pos[j],num_to_base_map[batch_ref[j]],num_to_base_map[pred1],num_to_base_map[pred2],min(999,-100*np.log10(1e-10+ 1-batch_probs[j,pred2])),'PASS','1/2', batch_dp[j], batch_freq[j])
                                f.write(s)

                    elif batch_pred_GT[j]==1 and batch_ref[j]!=batch_pred[j,-1] and batch_probs[j,batch_pred[j,-1]]>=0.5:
                        pred1=batch_pred[j,-1]
                        s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t.\tGT:DP:FQ\t%s:%d:%.4f\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred1], min(999,-100*np.log10(1e-10+ 1-batch_probs[j,pred1])),'PASS','1/1', batch_dp[j], batch_freq[j])
                        f.write(s)


                    neg_file.write('%d,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f\n' %(batch_pos[j],num_to_base_map[batch_ref[j]], batch_prob_GT[j,0], batch_probs[j,0], batch_probs[j,1], batch_probs[j,2], batch_probs[j,3], batch_dp[j], batch_freq[j]))
                    
                    
            print('%s: (%d/%d) regions completed.' %(str(datetime.datetime.now()), completed, total_regions),flush=True)
            f.flush()
            os.fsync(f.fileno())
            neg_file.flush()
            os.fsync(neg_file.fileno())
    
    output_file=os.path.join(vcf_path,'%s.snps' %prefix)
    run_cmd("bcftools sort %s.vcf|bgziptabix %s.vcf.gz" %(output_file, output_file))  
    
    return output_file
                    
