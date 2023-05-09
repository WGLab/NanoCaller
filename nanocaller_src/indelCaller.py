import time, os, datetime, queue
from tqdm import tqdm
import numpy as np
import multiprocessing as mp
import tensorflow as tf
from .model_architect_indel import *
from .generate_indel_pileups import get_indel_testing_candidates
from .utils import *
from multiprocessing import current_process


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


def indel_run(params, indel_dict, job_Q, counter_Q, indel_files_list):
    cur_p = current_process()
    curr_vcf_path=os.path.join(params['intermediate_indel_files_dir'],'%s.%d.indel.vcf' %(params['prefix'], cur_p._identity[0]))
    indel_files_list.append(curr_vcf_path)
    
    model_path=get_indel_model(params['indel_model'])
    if model_path==None:
        print('Invalid indel model name or path', flush=True)
        sys.exit(1)
    
    indel_model=Indel_model()    
    indel_model.load_weights(model_path).expect_partial()
    
    batch_size=100

    with open(curr_vcf_path,'w') as f:
        while len(indel_dict)>0 or not job_Q.empty():
            try:
                job=job_Q.get(block=False)
                chunk=job[1]
                chrom=chunk['chrom']
                pos, x0_test, x1_test, x2_test, alleles_seq, phase = get_indel_testing_candidates(params, chunk)
                
                prev=0
                prev_len=0

                if len(pos)>0:
                    for batch in range(int(np.ceil(len(x0_test)/batch_size))):
                        batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]

                        batch_x0 = x0_test[batch*batch_size:min((batch+1)*batch_size,len(x0_test))]
                        batch_x1 = x1_test[batch*batch_size:min((batch+1)*batch_size,len(x1_test))]
                        batch_x2 = x2_test[batch*batch_size:min((batch+1)*batch_size,len(x2_test))]
                        batch_alleles_seq = alleles_seq[batch*batch_size:min((batch+1)*batch_size,len(alleles_seq))]
                        batch_phase = phase[batch*batch_size:min((batch+1)*batch_size,len(phase))]

                        batch_x_all=np.hstack([batch_x0, batch_x1, batch_x2])

                        batch_prob_all= indel_model(batch_x_all)

                        batch_pred_all=np.argmax(batch_prob_all,axis=1)

                        qual_all=-10*np.log10(1e-6+1-batch_prob_all)

                        for j in range(len(batch_pred_all)):
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
                                                if batch_phase[j]:
                                                    s='%s\t%d\t.\t%s\t%s,%s\t%.2f\tPASS\t.\tGT:GQ:PS\t1|2:%.2f:%d\n' %(chrom, batch_pos[j], ref, alt1, alt2, q, gq, batch_phase[j])
                                                else:
                                                    s='%s\t%d\t.\t%s\t%s,%s\t%.2f\tPASS\t.\tGT:GQ\t1|2:%.2f\n' %(chrom, batch_pos[j], ref, alt1, alt2, q, gq)
                                                f.write(s)
                                                prev=batch_pos[j]+max(len(ref), len(alt1),len(alt2))

                                        elif allele0_data[0]:
                                            gq=-100*np.log10(1+1e-6-batch_prob_all[j,2])
                                            if batch_phase[j]:
                                                s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ:PS\t0|1:%.2f:%d\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], q, gq, batch_phase[j])
                                            else:
                                                s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ\t0|1:%.2f\n' %(chrom, batch_pos[j], allele0_data[0], allele0_data[1], q, gq)
                                            f.write(s)
                                            prev=batch_pos[j]+max(len(allele0_data[0]), len(allele0_data[1]))

                                        elif allele1_data[0]:
                                            gq=-100*np.log10(1+1e-6-batch_prob_all[j,2])
                                            if batch_phase[j]:
                                                s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ:PS\t1|0:%.2f:%d\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1], q, gq, batch_phase[j])
                                            else:
                                                s='%s\t%d\t.\t%s\t%s\t%.2f\tPASS\t.\tGT:GQ\t1|0:%.2f\n' %(chrom, batch_pos[j], allele1_data[0], allele1_data[1], q, gq)
                                            f.write(s)
                                            prev=batch_pos[j]+max(len(allele1_data[0]), len(allele1_data[1]))

                        batch_x0,batch_x1,batch_x2,batch_x=None,None,None,None
                    pos, x0_test, x1_test, x2_test=None,None,None,None

                f.flush()
                os.fsync(f.fileno())
                counter_Q.put(1)
            
            except queue.Empty:
                continue
    
    
def phase_run(contig_dict, params, indel_dict, job_Q, counter_Q, phased_snp_files_list):
    contig = contig_dict['name']
    start, end = contig_dict['start'], contig_dict['end']
    phase_dir = params['intermediate_phase_files_dir']
    input_snp_vcf = params['snp_vcf']
    input_contig_vcf = os.path.join(phase_dir, '%s.snps.unphased.vcf' %contig)
    raw_output_contig_vcf = os.path.join(phase_dir, '%s.snps.phased.raw.vcf' %contig)
    output_contig_vcf = os.path.join(phase_dir, '%s.snps.phased.vcf.gz' %contig)
    phased_bam = os.path.join(phase_dir, '%s.phased.bam' %contig)
    
    phased_snp_files_list.append(output_contig_vcf)
    
    enable_whatshap = '--distrust-genotypes --include-homozygous' if params['enable_whatshap'] else ''
    
    # extract region from VCF
    run_cmd('bcftools view %s -r %s -o %s' %(input_snp_vcf, contig, input_contig_vcf))
    
    #phase VCF
    run_cmd("whatshap phase %s %s -o %s -r %s --ignore-read-groups --chromosome %s %s" %(input_contig_vcf, params['sam_path'], raw_output_contig_vcf, params['fasta_path'], contig , enable_whatshap))


    run_cmd("bcftools view -e  'GT=\"0\\0\"' %s|bgziptabix %s" %(raw_output_contig_vcf, output_contig_vcf))

    run_cmd("whatshap haplotag --ignore-read-groups --ignore-linked-read --reference %s %s %s --regions %s:%d-%d --tag-supplementary -o - | samtools view -b -1 --write-index -o %s" %(params['fasta_path'], output_contig_vcf, params['sam_path'], contig, start, end, phased_bam), verbose=True)
    run_cmd("touch -c %s.csi" %phased_bam) 
    if params['mode']=='snps':
        counter_Q.put(1)
    else:
        indel_jobs=indel_dict[contig]
        for chunk in indel_jobs:
            chunk['sam_path']=phased_bam

            job_Q.put(('indel', chunk))

        indel_dict.pop(contig)
    
def caller(params, job_Q, counter_Q, indel_dict, phased_snp_files_list, indel_files_list):
    while len(indel_dict)>0 or not job_Q.empty():
        try:
            job=job_Q.get(block=False)
            if job[0]=='phase':
                phase_run(job[1], params, indel_dict, job_Q, counter_Q, phased_snp_files_list)

            elif job[0]=='indel':
                job_Q.put(job)
                indel_run(params, indel_dict, job_Q, counter_Q, indel_files_list)
                
        except queue.Empty:
            continue

def progress_bar(counter_Q, total_indel_jobs, mode, suppress_progress):
    if not suppress_progress:
        pbar = tqdm(total = total_indel_jobs)
        if mode=='snps':
            pbar.set_description("SNPs and BAM Phasing Progress")

        else:
            pbar.set_description("Indel Calling Progress")

        for item in iter(counter_Q.get, None):
            pbar.update()
        
def call_manager(params):
    pmanager = mp.Manager()
    
    indel_dict=pmanager.dict()
    job_Q=pmanager.Queue()
    counter_Q=pmanager.Queue()
    phased_snp_files_list=pmanager.list()
    indel_files_list=pmanager.list()
    
    contigs_list={}
    for x in params['regions_list']:
        if x[0] not in contigs_list:
            contigs_list[x[0]]={'name':x[0], 'start':x[1], 'end':x[2]}
        else:
            contigs_list[x[0]]['start']=min(x[1], contigs_list[x[0]]['start'])
            contigs_list[x[0]]['end']=max(x[2], contigs_list[x[0]]['end'])
    
    if params['mode']=='indels' or params['mode']=='all':
        params['intermediate_indel_files_dir']=os.path.join(params['vcf_path'], 'intermediate_indel_files')
        make_and_remove_path(params['intermediate_indel_files_dir'])
    
    if params['mode']=='snps' or params['mode']=='all':
        params['intermediate_phase_files_dir']=os.path.join(params['vcf_path'], 'intermediate_phase_files')
        make_and_remove_path(params['intermediate_phase_files_dir'])
        

    total_phase_jobs=0
    total_indel_jobs=0
    
    
    if params['mode']=='indels':
        for chunk in params['chunks_list']:
            chunk['sam_path']=params['sam_path']
            job_Q.put(('indel', chunk))
            total_indel_jobs+=1
    
    else:
        for contig in contigs_list:
            job_Q.put(('phase', contigs_list[contig]))
            total_phase_jobs+=1
            
        if params['mode']=='all':
            for chunk in params['chunks_list']:
                if chunk['chrom'] not in indel_dict:
                    indel_dict[chunk['chrom']]=pmanager.list()   
                indel_dict[chunk['chrom']].append(chunk)
                total_indel_jobs+=1
    
    total_jobs=total_phase_jobs if params['mode']=='snps' else total_indel_jobs
    
    tqdm_proc = mp.Process(target=progress_bar, args=(counter_Q, total_jobs, params['mode'], params['suppress_progress']))
    tqdm_proc.start()
    
    shared_var = (params, job_Q, counter_Q, indel_dict, phased_snp_files_list, indel_files_list)
    
    handlers = []
    for hid in range(params['cpu']):
        p = mp.Process(target=caller, args=shared_var);
        p.start();
        handlers.append(p);

        
    for job in handlers:
        job.join()
    
    counter_Q.put(None)
    tqdm_proc.join()
        
    output_files={'snps':None, 'indels':None, 'final':None}
    
    if params['mode']=='snps' or params['mode']=='all':
        output_files['snps']=os.path.join(params['vcf_path'],'%s.snps.phased.vcf.gz' %params['prefix'])
        
        phased_snps_files='\n'.join(phased_snp_files_list)
        cmd='bcftools concat -f - |bgziptabix %s' %(output_files['snps'])
        stream=Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = stream.communicate(input=phased_snps_files.encode())
        
    if params['mode']=='indels' or params['mode']=='all':
        raw_indel_vcf=os.path.join(params['intermediate_indel_files_dir'],'%s.raw.indel.vcf' %params['prefix'])
        output_files['indels']=os.path.join(params['vcf_path'],'%s.indels.vcf.gz' %params['prefix'])
        
        with open(raw_indel_vcf,'wb') as outfile:
            outfile.write(b'##fileformat=VCFv4.2\n')
            outfile.write(b'##FILTER=<ID=PASS,Description="All filters passed">\n')

            #loop through contigs
            for chrom in contigs_list:
                outfile.write(b'##contig=<ID=%s>\n' %bytes(chrom.encode('utf-8')))

            outfile.write(b'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            outfile.write(b'##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Probability">\n')
            outfile.write(b'##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">\n')
            outfile.write(b'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' %bytes(params['sample'].encode('utf-8')))

            for int_file in indel_files_list:
                with open(int_file,'rb') as fd:
                    shutil.copyfileobj(fd, outfile)
        
        print('\n%s: Compressing and indexing indel calls.' %str(datetime.datetime.now()))
        
        run_cmd(' bcftools sort %s|rtg RTG_MEM=2G vcfdecompose -i - -o - |rtg RTG_MEM=2G vcffilter -i - --non-snps-only -o  %s' %(raw_indel_vcf, output_files['indels']))
        
    if params['mode']=='all':
        final_vcf=os.path.join(params['vcf_path'],'%s.vcf.gz' %params['prefix'])
        output_files['final']=final_vcf
        
        run_cmd('bcftools concat %s %s -a | bgziptabix %s' %(output_files['snps'], output_files['indels'], output_files['final']), error=True)
        
        
    return output_files
