import time,os,copy,argparse,subprocess, sys, pysam, datetime
from tqdm import tqdm
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from .model_architect import *
from .model_architect_SNP_haploid import haploid_SNP_model
from .generate_SNP_pileups import get_snp_testing_candidates
from intervaltree import Interval, IntervalTree
from .utils import *
from multiprocessing import current_process

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
num_to_base_map={0:'A',1:'G',2:'T',3:'C'}

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
                'CLR-HG002':'release_data/clr_models/SNPs/HG002_giab-4.2.1/model-100',
                'haploid':'release_data/haploid_models/SNPs/CHM13/model.24-0.9985.h5'
               }

def get_SNP_model(snp_model):
    if snp_model in snp_model_dict:
        dirname = os.path.dirname(__file__)
        snp_model_path=os.path.join(dirname, snp_model_dict[snp_model])
        
    elif os.path.exists(snp_model):
        if os.path.isdir(snp_model):
            snp_model_path=glob.glob(os.path.join(snp_model,'*.index'))[0].rstrip('.index')
            
    else:
        return None,None
    
    coverage_path='%s.coverage' %snp_model_path
    
    if os.path.exists(coverage_path):
        train_coverage=float(open(coverage_path,'r').readlines()[0].rstrip('\n'))
    else:
        train_coverage=0
        
    return snp_model_path, train_coverage
    
def caller(params, chunks_Q, counter_Q, snp_files):
    cur_p = current_process()
    curr_vcf_path=os.path.join(params['intermediate_snp_files_dir'],'%s.%d.snps.vcf' %(params['prefix'], cur_p._identity[0]))
    snp_files.append(curr_vcf_path)
    
    n_input=[5,41,5]
    
    model_path, train_coverage=get_SNP_model(params['snp_model'])
    
    if model_path==None:
        print('Invalid SNP model name or path', flush=True)
        sys.exit(1)
            
    snp_model=SNP_model()    
    snp_model.load_weights(model_path).expect_partial()
    
    hap_train_coverage=30
    hap_model_path,_=get_SNP_model('haploid')
    hap_snp_model=haploid_SNP_model()    
    null_x=[np.zeros(5*41*5).reshape(1,5,41,5),np.zeros(4).reshape(1,4)]
    hap_snp_model(null_x);
    hap_snp_model.load_weights(hap_model_path)
    
    batch_size=1000
    
    with open(curr_vcf_path,'w') as f:
        while not chunks_Q.empty():
            chunk = chunks_Q.get(block=False)
            chrom=chunk['chrom']
            pos, test_ref, x_test, dp, freq, coverage, fwd_dp, rev_dp = get_snp_testing_candidates(params, chunk)
            
            if len(pos)!=0:                
                if chunk['ploidy']=='diploid':
                    test_ref=test_ref.astype(np.float16)
                    x_test=x_test.astype(np.float32)
                    
                    if params['disable_coverage_normalization']:
                        x_test[:,1:,:,:4]=x_test[:,1:,:,:4]*(train_coverage/dp[:,np.newaxis,np.newaxis,np.newaxis])
                    else:
                        x_test[:,1:,:,:4]=x_test[:,1:,:,:4]*(train_coverage/coverage)

                    for batch in range(int(np.ceil(len(x_test)/batch_size))):
                        batch_freq=freq[batch*batch_size:min((batch+1)*batch_size,len(freq))]
                        batch_dp=dp[batch*batch_size:min((batch+1)*batch_size,len(dp))]
                        batch_fwd_dp=fwd_dp[batch*batch_size:min((batch+1)*batch_size,len(fwd_dp))]
                        batch_rev_dp=rev_dp[batch*batch_size:min((batch+1)*batch_size,len(rev_dp))]
                        
                        batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                        batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]

                        batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]

                        batch_coverage=np.mean(batch_dp)

                        batch_prob_A, batch_prob_G, batch_prob_T, batch_prob_C, batch_prob_GT = snp_model([batch_x, batch_ref[:,0][:,np.newaxis], batch_ref[:,1][:,np.newaxis], batch_ref[:,2][:,np.newaxis], batch_ref[:,3][:,np.newaxis]])

                        batch_pred_GT=np.argmax(batch_prob_GT,axis=1)

                        batch_probs=np.hstack([batch_prob_A[:,1][:,np.newaxis], batch_prob_G[:,1][:,np.newaxis], batch_prob_T[:,1][:,np.newaxis], batch_prob_C[:,1][:,np.newaxis]])


                        batch_pred=np.argsort(batch_probs,axis=1)

                        batch_ref=np.argmax(batch_ref,1)

                        batch_pred_GT=np.sum(batch_probs>=0.5,axis=1)

                        sort_probs=np.sort(batch_probs,axis=1)

                        for j in range(len(batch_pred_GT)):
                            info_field='PR='+','.join("{:.4f}".format(x) for x in batch_probs[j,[0,3,1,2]])+ ";FQ={:.4f}".format(batch_freq[j])
                            ref_dp=(batch_fwd_dp[j][batch_ref[j]], batch_rev_dp[j][batch_ref[j]])
                            
                            if batch_pred_GT[j]>=2: # if het
                                    pred1,pred2=batch_pred[j,-1],batch_pred[j,-2]
                                    if pred1==batch_ref[j]:
                                        alt_dp=(batch_fwd_dp[j][pred2], batch_rev_dp[j][pred2])
                                        
                                        s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:%.4f:%d,%d:%d,%d:%d,%d\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred2], min(99,-10*np.log10(1e-10+ 1-batch_probs[j,pred2])),'PASS',info_field, '0/1', batch_dp[j], sum(alt_dp)/batch_dp[j], sum(ref_dp), sum(alt_dp), ref_dp[0], alt_dp[0], ref_dp[1], alt_dp[1])
                                        f.write(s)

                                    elif pred2==batch_ref[j] and batch_probs[j,pred2]>=0.5:
                                        alt_dp=(batch_fwd_dp[j][pred1], batch_rev_dp[j][pred1])
                                        s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:%.4f:%d,%d:%d,%d:%d,%d\n' %(chrom,batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred1], min(99,-10*np.log10(1e-10+ 1-batch_probs[j,pred2])),'PASS',info_field, '0/1', batch_dp[j], sum(alt_dp)/batch_dp[j], sum(ref_dp), sum(alt_dp), ref_dp[0], alt_dp[0], ref_dp[1], alt_dp[1])
                                        f.write(s)

                                    elif pred2!=batch_ref[j] and pred1!=batch_ref[j] and batch_probs[j,pred2]>=0.5:
                                        alt1_dp=(batch_fwd_dp[j][pred1], batch_rev_dp[j][pred1])
                                        alt2_dp=(batch_fwd_dp[j][pred2], batch_rev_dp[j][pred2])
                                        s='%s\t%d\t.\t%s\t%s,%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:%.4f,%.4f:%d,%d,%d:%d,%d,%d:%d,%d,%d\n' %\
                            (chrom,batch_pos[j],num_to_base_map[batch_ref[j]],num_to_base_map[pred1],num_to_base_map[pred2],min(99,-10*np.log10(1e-10+ 1-batch_probs[j,pred2])),'PASS',info_field,'1/2', batch_dp[j], sum(alt1_dp)/batch_dp[j],sum(alt2_dp)/batch_dp[j], sum(ref_dp), sum(alt1_dp), sum(alt2_dp), ref_dp[0], alt1_dp[0], alt2_dp[0], ref_dp[1], alt1_dp[1], alt2_dp[1])
                                        f.write(s)

                            elif batch_pred_GT[j]==1 and batch_ref[j]!=batch_pred[j,-1] and batch_probs[j,batch_pred[j,-1]]>=0.5:
                                pred1=batch_pred[j,-1]
                                alt_dp=(batch_fwd_dp[j][pred1], batch_rev_dp[j][pred1])
                                s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:%.4f:%d,%d:%d,%d:%d,%d\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred1], min(99,-10*np.log10(1e-10+ 1-batch_probs[j,pred1])),'PASS',info_field,'1/1', batch_dp[j], sum(alt_dp)/batch_dp[j], sum(ref_dp), sum(alt_dp), ref_dp[0], alt_dp[0], ref_dp[1], alt_dp[1])
                                f.write(s)
                                
                            else:
                                if batch_pred_GT[j]==1 and batch_ref[j]==batch_pred[j,-1]:
                                    pred1=batch_pred[j,-1]
                                    s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:.:.:.:.\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], '.', min(99,-10*np.log10(1e-10+ 1-batch_probs[j,pred1])),'REF',info_field,'./.', batch_dp[j])
                                    f.write(s)
                                else:
                                    s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:.:.:.:.\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], '.', 0,'LOW',info_field,'./.', batch_dp[j])
                                    f.write(s)
                                    

                elif chunk['ploidy']=='haploid':
                    test_ref=test_ref.astype(np.float16)
                    x_test=x_test.astype(np.float32)
                    
                    if params['disable_coverage_normalization']:
                        x_test[:,1:,:,:4]=x_test[:,1:,:,:4]*(hap_train_coverage/dp[:,np.newaxis,np.newaxis,np.newaxis])
                    else:
                        x_test[:,1:,:,:4]=x_test[:,1:,:,:4]*(hap_train_coverage/coverage)
                        
                    for batch in range(int(np.ceil(len(x_test)/batch_size))):
                        batch_freq=freq[batch*batch_size:min((batch+1)*batch_size,len(freq))]
                        batch_dp=dp[batch*batch_size:min((batch+1)*batch_size,len(dp))]
                        batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                        batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]

                        batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]    

                        batch_probs=hap_snp_model([batch_x,batch_ref]).numpy()
                        batch_ref=np.argmax(batch_ref,1)
                        batch_pred=np.argmax(batch_probs,1)

                        for j in range(len(batch_pred)):
                            pred=batch_pred[j]
                            
                            if pred!=batch_ref[j]:
                                info_field='PR='+','.join("{:.4f}".format(x) for x in batch_probs[j,[0,3,1,2]])+ ";FQ={:.4f}".format(batch_freq[j])
                                s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:%.4f:.:.:.\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred], min(999,-100*np.log10(1e-10+ 1-batch_probs[j,pred])),'PASS',info_field,'1/1', batch_dp[j], batch_freq[j])
                                f.write(s)
                                
                            else:
                                info_field='PR='+','.join("{:.4f}".format(x) for x in batch_probs[j,[0,3,1,2]])+ ";FQ={:.4f}".format(batch_freq[j])
                                s='%s\t%d\t.\t%s\t%s\t%.3f\t%s\t%s\tGT:DP:VF:AD:ADF:ADR\t%s:%d:%.4f:.:.:.\n' %(chrom, batch_pos[j], num_to_base_map[batch_ref[j]], num_to_base_map[pred], min(999,-100*np.log10(1e-10+ 1-batch_probs[j,pred])),'REF',info_field,'1/1', batch_dp[j], batch_freq[j])
                                f.write(s)

                f.flush()

                os.fsync(f.fileno())
            counter_Q.put(1)

def progress_bar(counter_Q, total_snp_jobs, suppress_progress):
    if not suppress_progress:
        pbar = tqdm(total = total_snp_jobs)
        pbar.set_description("SNP Calling Progress")
        for item in iter(counter_Q.get, None):
            pbar.update()

        
def call_manager(params):
    
    pmanager = mp.Manager()
    chunks_Q = pmanager.Queue()
    snp_files = pmanager.list()    
    counter_Q = pmanager.Queue()
    
    chrom_list=set(x[0] for x in params['regions_list'])
    
    total_snp_jobs=0
    for chunk in params['chunks_list']:
        chunks_Q.put(chunk)
        total_snp_jobs+=1
        
    
    params['intermediate_snp_files_dir']=os.path.join(params['vcf_path'], 'intermediate_snp_files')
    make_and_remove_path(params['intermediate_snp_files_dir'])
    
    
    tqdm_proc = mp.Process(target=progress_bar, args=(counter_Q, total_snp_jobs, params['suppress_progress']))
    tqdm_proc.start()
    
    
    shared_var = (params, chunks_Q, counter_Q, snp_files)
    handlers = []
    for hid in range(params['cpu']):
        p = mp.Process(target=caller, args=shared_var);
        p.start();
        handlers.append(p);

        
    for job in handlers:
        job.join()
    
    
    counter_Q.put(None)
    tqdm_proc.join()
    
    all_snps_file_path=os.path.join(params['vcf_path'],'%s.unfiltered.snps.vcf.gz' %params['prefix'])
    pass_snps_file_path=os.path.join(params['vcf_path'],'%s.snps.vcf.gz' %params['prefix'])
    
    temp_output_file_path=os.path.join(params['intermediate_snp_files_dir'],'combined.snps.vcf')
    
    print('\n%s: Combining SNP calls.' %str(datetime.datetime.now()))
    
    with open(temp_output_file_path,'wb') as outfile:
        outfile.write(b'##fileformat=VCFv4.2\n')
        outfile.write(b'##FILTER=<ID=PASS,Description="All filters passed">\n')
        outfile.write(b'##FILTER=<ID=LOW,Description="All alleles have probability less than 50%.">\n')
        outfile.write(b'##FILTER=<ID=REF,Description="Homozygous Reference. Only reference allele has greater than 50% probability. All alternative alleles having probability less than 50%.">\n')
        #loop through contigs
        for chrom in chrom_list:
            outfile.write(b'##contig=<ID=%s>\n' %bytes(chrom.encode('utf-8')))

        outfile.write(b'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        outfile.write(b'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        outfile.write(b'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
        outfile.write(b'##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Allelic depths on forward strand for the ref and alt alleles in the order listed">\n')
        outfile.write(b'##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Allelic depths on reverse strand for the ref and alt alleles in the order listed">\n')
        outfile.write(b'##FORMAT=<ID=VF,Number=A,Type=Float,Description="Alternative allele frequency in the order listed">\n')
        outfile.write(b'##INFO=<ID=PR,Number=4,Type=Float,Description="Probability of presence of alleles A, C, G and T, in the given order. Probability of each base is out of 1, independent of each other.">\n')
        outfile.write(b'##INFO=<ID=FQ,Number=1,Type=Float,Description="Maximum frequency of non-reference base.">\n')
        
        outfile.write(b'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' %bytes(params['sample'].encode('utf-8')))

        for int_file in snp_files:
            with open(int_file,'rb') as fd:
                shutil.copyfileobj(fd, outfile)
    
    print('\n%s: Compressing and indexing SNP calls.' %str(datetime.datetime.now()))
    
    run_cmd("bcftools sort %s --temp-dir %s/|bgzip > %a && tabix -fp vcf --csi %a" %(temp_output_file_path, params['intermediate_snp_files_dir'], all_snps_file_path, all_snps_file_path), error=True)
    run_cmd("bcftools view %s -f PASS|bgzip > %a && tabix -fp vcf --csi %a" %(all_snps_file_path, pass_snps_file_path, pass_snps_file_path), error=True)
    
    return pass_snps_file_path
    
