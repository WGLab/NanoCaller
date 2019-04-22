import pysam, sys
import pandas as pd
import numpy as np
import multiprocessing as mp
import zipfile
import io,time,os

def get_candidate(chr_num,start,end,dl_lim,threshold,samfile,fastafile):
    # this function finds all the candidate variant sites in the range given by chr_num:start-end
    output=[]
    for pcol in samfile.pileup(chr_num,start-1,end,min_base_quality=0,stepper='nofilter',\
                                       flag_filter=0x800,truncate=True):
        
        #if number of reads<dl_lim or reference base percentage is less than threshold discard the site                
        n=pcol.get_num_aligned()
        xx=[x.upper() for x in pcol.get_query_sequences()]
        r=fastafile.fetch('chr20',start=pcol.pos,end=pcol.pos+1)
        #if number of reads<dl_lim or reference base percentage is less than threshold discard the site
        if n>=dl_lim and xx.count(r)/n<=threshold:
            output.append(pcol.pos+1)
    return output

def create_pileup(options):
    #this function creates pileup of all the reference positions in the range chr_num:start-end
    chr_num,start,end,sam_path,fasta_path=options
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)
    
    mapping={'*':4,'A':0,'G':1,'T':2,'C':3,'N':5}

    d={x:{} for x in range(start,end)}
    features={x:{} for x in range(start,end)}
    rlist=[mapping[s] for s in fastafile.fetch('chr20',start-1,end-1)]
    for pcol in samfile.pileup(chr_num,start-1,end-1,min_base_quality=0,stepper='nofilter',\
                                       flag_filter=0x800,truncate=True):
        name=pcol.get_query_names()
        seq=pcol.get_query_sequences()
        qual=pcol.get_query_qualities()

        d[pcol.pos+1]={n:s.upper() if len(s)>0 else '*' for (n,s) in zip(name,seq)}
        features[pcol.pos+1]={n:q for (n,q) in zip(name,qual)}
    
    #create pileup dataframe
    p_df=pd.DataFrame.from_dict(d)
    p_df.fillna('N',inplace=True)
    p_df=p_df.applymap(lambda x: mapping[x])
    p_mat=np.array(p_df)
    
    
    #generatre quality dataframe
    f_df=pd.DataFrame.from_dict(features)
    f_df.fillna(0,inplace=True)
    f_df=f_df.reindex(p_df.index)
    
    #create reference match matrix containing +1 for mathc and -1 for non-match
    ref_match=(p_mat==rlist)
    tmp2=2*np.array(ref_match)[:,:,np.newaxis]-1
    
    #encode bases as one hot vectors, and multiply with +1 if it matches with reference and -1 otherwise
    #combine all matrices
    tmp=np.dstack([(p_mat==i) for i in range(4)])
    data=np.multiply(tmp,tmp2).astype(np.int8)
    data=np.dstack((data,np.array(f_df).astype(np.int8)))
    
    
    return (data,list(p_df.index),str(start)+'-'+str(end))

def generate(chr_num,start,end,sam_path,fasta_path,out_path):
    #this function calls both get_candidates and create_pileup with given parameters
    du_lim=32
    dl_lim=12
    threshold=0.7
    samfile = pysam.Samfile(sam_path, "rb")
    fastafile=pysam.FastaFile(fasta_path)
    
    #get candidates and save them
    candidates=get_candidate(chr_num,start,end,dl_lim,threshold,samfile,fastafile)
    name=chr_num+'_'+str(start)+'_'+str(end)+'_'+'candidates'
    np.savez_compressed(os.path.join(out_path,name), candidates=candidates)
    
    #generate pileups by breaking the range into 1mb regions. Each 1mb region is further split in ten and
    #these pileups are generated in parallel. Each 1mb region is stored in one file.
    print('starting pileups')
    for mbase in range(start,end,int(1e6)):
        print('starting pool:'+str(mbase))
        name=chr_num +'_'+ str(mbase) +'_'+'.npz'
        pool = mp.Pool(processes=mp.cpu_count())
        in_dict=[[chr_num,k-100,k+100100,sam_path,fasta_path] for k in range(mbase,mbase+int(1e6),100000)]
        results = pool.map(create_pileup, in_dict)
        pileups={}
        for res in results:
            pileups[res[2]]=res[0]
            pileups[res[2]+'-rnames']=res[1]
        saveCompressed(open(os.path.join(out_path,name), 'wb'),pileups)
        print('finishing pool:'+str(mbase))
    
    
def saveCompressed(fh, namedict):
     with zipfile.ZipFile(fh, mode="w", compression=zipfile.ZIP_DEFLATED,
                          allowZip64=True) as zf:
        for k, v in namedict.items():
             with zf.open(k + '.npy', 'w', force_zip64=True) as buf:
                    np.lib.npyio.format.write_array(buf,np.asanyarray(v),allow_pickle=False)
    

    
if __name__ == '__main__':
    if len(sys.argv)!=7:
        print('Wrong number of inputs')
        sys.exit(0)
    
    chr_num,start,end,sam_path,fasta_path,out_path=sys.argv[1:7]
    generate(chr_num,int(start),int(end),sam_path,fasta_path,out_path)
    
    
