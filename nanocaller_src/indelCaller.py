import time,os,copy,argparse,subprocess, datetime

import numpy as np
import multiprocessing as mp
import tensorflow as tf
from .model_architect_indel import *
from Bio import pairwise2
from .generate_indel_pileups import get_indel_testing_candidates
from .utils import *

rev_gt_map={0:'hom-ref', 1:'hom-alt', 2:'het-ref', 3:'het-alt'}
rev_base_map={0:'A',1:'G',2:'T',3:'C',4:'-'}

indel_model_dict={'NanoCaller1':'release_data/ONT_models/indels/NanoCaller1_beta/model-30',
                  'NanoCaller3':'release_data/hifi_models/indels/NanoCaller3_beta/model-25',
                  'ONT-HG001':'release_data/ONT_models/indels/HG001_guppy4.2_giab-3.3.2/model-100',
                  'ONT-HG002':'release_data/ONT_models/indels/HG002_guppy4.2_giab-4.2.1/model-100',
                  'CCS-HG001':'release_data/hifi_models/indels/HG001_giab-3.3.2/model-100',
                  'CCS-HG002':'release_data/hifi_models/indels/HG002_giab-4.2.1/model-100'}

def get_indel_model(indel_model):
    if os.path.exists(indel_model):
        if os.path.isdir(indel_model):
            indel_model_path=glob.glob(os.path.join(indel_model,'*.index'))[0].rstrip('.index')
    
    elif indel_model in indel_model_dict:
        dirname = os.path.dirname(__file__)
        indel_model_path=os.path.join(dirname, indel_model_dict[indel_model])
        
    else:
        return None

    return indel_model_path


def test_model(params,pool):
    print('%s: Indel calling started.' %(str(datetime.datetime.now())), flush=True)
    
    rev_allele_map={0:'N',1:'D',2:'I'}
    
    vcf_path,prefix=params['vcf_path'], params['prefix']
    chrom,start,end=params['chrom'], params['start'],params['end']
    
    model_path=get_indel_model(params['indel_model'])
    
    if model_path==None:
        print('Invalid indel model name or path', flush=True)
        return
    
    indel_model=Indel_model()    
    indel_model.load_weights(model_path)
    
    batch_size=100

    neg_file=open(os.path.join(vcf_path,'%s.indel_stats' %prefix),'w')
    neg_file.write('pos,ref\n')
    
    
    outfile=os.path.join(vcf_path,'%s.indels' %prefix)
    
    with open('%s.raw.vcf' %outfile,'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Probability">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s\n' %params['sample'])
        
        in_dict_list=[]
        
        for mbase in range(start, end, 50000):
            d = copy.deepcopy(params)
            d['start']=mbase
            d['end']=min(end,mbase+50000)
            in_dict_list.append(d)
        
        result=pool.imap_unordered(get_indel_testing_candidates, in_dict_list)
        
        total_regions=len(in_dict_list)
        completed=0
        for res in result:            
            pos, x0_test, x1_test, x2_test, alleles_seq=res
            completed+=1
            
            prev=0
            prev_len=0
            
            if len(pos)==0:
                continue

            for batch in range(int(np.ceil(len(x0_test)/batch_size))):
                batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]

                batch_x0 = x0_test[batch*batch_size:min((batch+1)*batch_size,len(x0_test))]
                batch_x1 = x1_test[batch*batch_size:min((batch+1)*batch_size,len(x1_test))]
                batch_x2 = x2_test[batch*batch_size:min((batch+1)*batch_size,len(x2_test))]
                batch_alleles_seq = alleles_seq[batch*batch_size:min((batch+1)*batch_size,len(alleles_seq))]

                batch_x_all=np.hstack([batch_x0, batch_x1, batch_x2])
                
                batch_prob_all= indel_model(batch_x_all)
                
                batch_pred_all=np.argmax(batch_prob_all,axis=1)
                
                qual_all=-10*np.log10(1e-6+1-batch_prob_all)
                
                for j in range(len(batch_pred_all)):
                    neg_file.write('%s,%d,%s,%.4f,%.4f,%.4f,%.4f\n' %(chrom,batch_pos[j], rev_gt_map[batch_pred_all[j]], batch_prob_all[j,0], batch_prob_all[j,1], batch_prob_all[j,2], batch_prob_all[j,3]))
                    
                    q=min(999,qual_all[j,batch_pred_all[j]])
                    if batch_pos[j]>prev:
                        
                        if batch_prob_all[j,0]<=0.95:
                            
                            q=-100*np.log10(1e-6+batch_prob_all[j,0])
                            allele0_data, allele1_data,allele_total_data= batch_alleles_seq[j]
                            
                            if batch_pred_all[j]==1 and allele_total_data[0]:
                                        gq=-100*np.log10(1+1e-6-batch_prob_all[j,1])
                                        s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ\t1/1:%.2f\n' %(chrom, batch_pos[j], allele_total_data[0], allele_total_data[1], q, gq)
                                        f.write(s)
                                        prev=batch_pos[j]+max(len(allele_total_data[0]), len(allele_total_data[1]))
                            
                            else:
                                if allele0_data[0] and allele1_data[0]:
                                    
                                    if allele0_data[0]==allele1_data[0] and allele0_data[1]==allele1_data[1]:
                                        gq=-100*np.log10(1+1e-6-batch_prob_all[j,1])
                                        s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ\t1/1:%.2f\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], q, gq)
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
                                        s='%s\t%d\t.\t%s\t%s,%s\t%.2f\tPASS\t.\tGT:GQ\t1|2:%.2f\n' %(chrom, batch_pos[j], ref, alt1, alt2, q, gq)
                                        f.write(s)
                                        prev=batch_pos[j]+max(len(ref), len(alt1),len(alt2))

                                elif allele0_data[0]:
                                    gq=-100*np.log10(1+1e-6-batch_prob_all[j,2])
                                    s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ\t0|1:%.2f\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], q, gq)
                                    f.write(s)
                                    prev=batch_pos[j]+max(len(allele0_data[0]), len(allele0_data[1]))

                                elif allele1_data[0]:
                                    gq=-100*np.log10(1+1e-6-batch_prob_all[j,2])
                                    s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ\t1|0:%.2f\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1], q, gq)
                                    f.write(s)
                                    prev=batch_pos[j]+max(len(allele1_data[0]), len(allele1_data[1]))

                batch_x0,batch_x1,batch_x2,batch_x=None,None,None,None
            pos, x0_test, x1_test, x2_test=None,None,None,None
            
            print('%s: (%d/%d) regions completed.' %(str(datetime.datetime.now()), completed, total_regions),flush=True)
            f.flush()
            os.fsync(f.fileno())
            neg_file.flush()
            os.fsync(neg_file.fileno())
            
            
    neg_file.close()
    
    run_cmd("bcftools sort %s.raw.vcf|bgziptabix %s.raw.vcf.gz" %(outfile, outfile))

    return outfile
