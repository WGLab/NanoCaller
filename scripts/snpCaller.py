import time,os,copy,argparse,subprocess, sys
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect import *
from generate_SNP_pileups import generate


config = tf.ConfigProto()
config.gpu_options.allow_growth = True

num_to_base_map={0:'A',1:'G',2:'T',3:'C'}
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

def test_model(params,pool):
    model,vcf_path,prefix= params['model'],params['vcf_path'], params['prefix']
    chrom,start,end=params['chrom'], params['start'],params['end']
    
    if params['supplementary']:
        flag=0x4|0x100|0x200|0x400
    else:
        flag=0x4|0x100|0x200|0x400|0x800
        
    stream = os.popen("samtools depth %s -r %s:%d-%d -G %d|awk '{if ($3>=%d){sum+=$3; cnt+=1}}END{if(cnt>0){print sum/cnt}else{print 0}}'" %(params['sam_path'], chrom,start,end,flag, params['mincov']))
    coverage=float(stream.read())
    print('Coverage=%.2fx' %coverage, flush=True)
    if coverage==0:
        print('No coverage found', flush=True)
        sys.exit(0)
    
    
    n_input=[5,41,5]
    cpu=params['cpu']
    tf.reset_default_graph()
    
    weights,biases,tensors=get_tensors(n_input,0.0)
    (x,GT_label,A_label, G_label, T_label, C_label,GT_score, A_score, G_score, T_score, C_score, accuracy_GT, accuracy_A,  accuracy_G,  accuracy_T,  accuracy_C, prediction_accuracy_GT, prediction_accuracy_A,  prediction_accuracy_G,  prediction_accuracy_T,  prediction_accuracy_C, prediction_GT, prediction_A,  prediction_G,  prediction_T,  prediction_C, accuracy, cost, optimizer,  cost_GT, cost_A, cost_G, cost_T, cost_C,A_ref,G_ref,T_ref,C_ref,prob_GT,prob_A,prob_G,prob_T,prob_C,keep)=tensors

    

    
    dirname = os.path.dirname(__file__)
    if model== 'NanoCaller1':
        model_path=os.path.join(dirname, 'release_data/NanoCaller1/model-rt-1')
        train_coverage=43
    
    elif model== 'NanoCaller2':
        model_path=os.path.join(dirname, 'release_data/NanoCaller2/model-rt-1')
        train_coverage=62
        
    elif model== 'NanoCaller3':
        model_path=os.path.join(dirname, 'release_data/NanoCaller3/model-rt-100')
        train_coverage=28
        
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
        
        for mbase in range(start,end,int(1e7)):
            d = copy.deepcopy(params)
            d['start']=mbase
            d['end']=min(end,mbase+int(1e7))
            pos,x_test,test_ref,dp,freq=generate(d,pool)
            
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
            f.flush()
            os.fsync(f.fileno())
            neg_file.flush()
            os.fsync(neg_file.fileno())
    
    output_file=os.path.join(vcf_path,'%s.snps' %prefix)
    stream=os.popen("bcftools sort %s.vcf|bgziptabix %s.vcf.gz" %(output_file, output_file))
    stream.read()
    
    
    
    return output_file
                    
