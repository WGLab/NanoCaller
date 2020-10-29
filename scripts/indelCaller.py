import time,os,copy,argparse,subprocess, datetime
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from model_architect_indel import *
from Bio import pairwise2
import generate_indel_pileups

config =  tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

rev_gt_map={0:'hom-ref', 1:'hom-alt', 2:'het-ref', 3:'het-alt'}
rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

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
    saver.restore(sess, os.path.join(dirname, 'release_data/NanoCaller_indels/ont/model-30' if params['seq']=='ont' else  'release_data/NanoCaller_indels/pacbio/model-25'))
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
        
        in_dict_list=[]
        
        for mbase in range(start, end, 50000):
            d = copy.deepcopy(params)
            d['start']=mbase
            d['end']=min(end,mbase+50000)
            in_dict_list.append(d)
        
        result=pool.imap_unordered(generate_indel_pileups.get_indel_testing_candidates, in_dict_list)
        
        for res in result:            
            pos, x0_test, x1_test, x2_test, alleles_seq=res
            if len(pos)==0:
                continue

            for batch in range(int(np.ceil(len(x0_test)/batch_size))):
                batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]

                batch_x0 = x0_test[batch*batch_size:min((batch+1)*batch_size,len(x0_test))]
                batch_x1 = x1_test[batch*batch_size:min((batch+1)*batch_size,len(x1_test))]
                batch_x2 = x2_test[batch*batch_size:min((batch+1)*batch_size,len(x2_test))]
                batch_alleles_seq = alleles_seq[batch*batch_size:min((batch+1)*batch_size,len(alleles_seq))]

                batch_x_all=np.hstack([batch_x0, batch_x1, batch_x2])
                
                batch_prob_all= sess.run(prob, feed_dict={x2: batch_x_all, rate:0.0})
                
                batch_pred_all=np.argmax(batch_prob_all,axis=1)
                
                qual_all=-10*np.log10(1e-6+1-batch_prob_all)
                
                for j in range(len(batch_pred_all)):
                    neg_file.write('%s,%d,%s,%.4f,%.4f,%.4f,%.4f\n' %(chrom,batch_pos[j], rev_gt_map[batch_pred_all[j]], batch_prob_all[j,0], batch_prob_all[j,1], batch_prob_all[j,2], batch_prob_all[j,3]))
                    
                    q=min(999,qual_all[j,batch_pred_all[j]])
                    if batch_pos[j]>prev:
                        
                        if batch_prob_all[j,0]<=0.95:
                            
                            q=-100*np.log10(1e-6+batch_prob_all[j,0])
                            allele0_data, allele1_data,allele_total_data= batch_alleles_seq[j]
                            
                            if batch_pred_all[j]==1:
                                    if allele_total_data[0]:
                                        gq=-100*np.log10(1+1e-6-batch_prob_all[j,1])
                                        s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP:GQ\t1/1:%s:%.2f\n' %(chrom, batch_pos[j], allele_total_data[0], allele_total_data[1], q, rev_gt_map[batch_pred_all[j]], gq)
                                        f.write(s)
                                        prev=batch_pos[j]+max(len(allele_total_data[0]), len(allele_total_data[1]))
                            
                            else:
                                if allele0_data[0] and allele1_data[0]:
                                    
                                    if allele0_data[0]==allele1_data[0] and allele0_data[1]==allele1_data[1]:
                                        gq=-100*np.log10(1+1e-6-batch_prob_all[j,1])
                                        s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP:GQ\t1/1:%s:%.2f\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1],q,rev_gt_map[batch_pred_all[j]], gq)
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
                                        gq=-100*np.log10(1+1e-6-batch_prob_all[j,3])
                                        s='%s\t%d\t.\t%s\t%s,%s\t%.2f\tPASS\t.\tGT:TP:GQ\t1|2:%s:%.2f\n' %(chrom, batch_pos[j], ref, alt1, alt2, q, rev_gt_map[batch_pred_all[j]], gq)
                                        f.write(s)
                                        prev=batch_pos[j]+max(len(ref), len(alt1),len(alt2))

                                elif allele0_data[0]:
                                    gq=-100*np.log10(1+1e-6-batch_prob_all[j,2])
                                    s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP:GQ\t0|1:%s:%.2f\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], q, rev_gt_map[batch_pred_all[j]], gq)
                                    f.write(s)
                                    prev=batch_pos[j]+max(len(allele0_data[0]), len(allele0_data[1]))

                                elif allele1_data[0]:
                                    gq=-100*np.log10(1+1e-6-batch_prob_all[j,2])
                                    s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:TP:GQ\t1|0:%s:%.2f\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1], q, rev_gt_map[batch_pred_all[j]], gq)
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