from subprocess import PIPE, Popen
import os, shutil

def run_cmd(cmd, verbose=False, output=False,error=False):
    stream=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = stream.communicate()
    
    stdout=stdout.decode('utf-8')
    stderr=stderr.decode('utf-8')
    
    if stderr:
        print(stderr, flush=True)
    
    if verbose:
        print(stdout, flush=True)
        
        
    if output:
        return stdout
    if error:
        return stderr
    

def remove_path(path):
    if os.path.exists(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)

def make_and_remove_path(path):
    remove_path(path)
    os.makedirs(path)
    

def get_coverage(params,pool):
    chrom,start,end=params['chrom'], params['start'],params['end']
    
    regions='-b %s' %params['include_bed'] if params['include_bed'] else ''
    
    if params['supplementary']:
        flag=0x4|0x100|0x200|0x400
    else:
        flag=0x4|0x100|0x200|0x400|0x800
        
    if params['cpu']==1:
        results=[get_chunk_coverage("""samtools depth %s -r %s:%d-%d -G %d %s|awk -v sum=0 -v count=0 '{if ($3>=%d){sum+=$3; count+=1}}END{print sum "," count}'""" %(params['sam_path'], chrom, start, end, flag, regions, params['mincov']))]

    else:
        coverage_regions=[]
    
        for mbase in range(start,end,200000):
            coverage_regions.append("""samtools depth %s -r %s:%d-%d -G %d %s|awk -v sum=0 -v count=0 '{if ($3>=%d){sum+=$3; count+=1}}END{print sum "," count}'""" %(params['sam_path'], chrom, mbase, min(end,mbase+200000), flag, regions, params['mincov']))
            
        results=pool.map(get_chunk_coverage, coverage_regions)
            
    total_bases=sum([x[0] for x in results])
    total_positions=sum([x[1] for x in results])

    
    coverage=total_bases/total_positions if total_bases*total_positions!=0 else 0
    return coverage

def get_chunk_coverage(cmd):
    stream=run_cmd(cmd, output= True).rstrip('\n').split(',')
    return (int(stream[0]),int(stream[1]))
