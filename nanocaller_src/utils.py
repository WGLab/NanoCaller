from subprocess import PIPE, Popen
import os, shutil, pysam

def get_regions_list(args):
    regions_list=[]
    if args.wgs_contigs:
        sam_file=pysam.Samfile(args.bam)
        
        for contig in list(range(1,23)) + ['X','Y']:
            contig='chr'+str(contig) if args.wgs_contigs=='chr1-22XY' else str(contig)
            if sam_file.is_valid_reference_name(contig):
                regions_list.append([contig, 1, sam_file.get_reference_length(contig), 'haploid' if args.haploid_genome else 'diploid'])
            
    elif args.regions:
        sam_file=pysam.Samfile(args.bam)
        for r in args.regions:
            r2=r.split(':')
            if len(r2)==1:
                r2=r2[0]
                if sam_file.is_valid_reference_name(r2):
                    regions_list.append([r2, 1, sam_file.get_reference_length(r2),'haploid' if args.haploid_genome else 'diploid'])
                else:
                    print('\n%s: Contig %s not present in the BAM file.'  %(str(datetime.datetime.now()), r2), flush=True)
            
            elif len(r2)==2:
                cord=r2[1].split('-')
                if len(cord)==2:
                    regions_list.append([r2[0], int(cord[0]), int(cord[1]), 'haploid' if args.haploid_genome else 'diploid'])
                else:
                    print('\n%s: Invalid region %s.'  %(str(datetime.datetime.now()), r), flush=True)
                    
            else:
                print('\n%s: Invalid region %s.'  %(str(datetime.datetime.now()), r), flush=True)
                
    elif args.bed:
        sam_file=pysam.Samfile(args.bam)
        with open(args.bed, 'r') as bed_file:
            for line in bed_file:
                line=line.rstrip('\n').split()
                if sam_file.is_valid_reference_name(line[0]):
                    regions_list.append([line[0], int(line[1]), int(line[2]), 'haploid' if args.haploid_genome else 'diploid'])
                else:
                    print('\n%s: Contig %s not present in the BAM file.'  %(str(datetime.datetime.now()), line[0]), flush=True)
    else:
        sam_file=pysam.Samfile(args.bam)
        regions_list=[[r, 1, sam_file.get_reference_length(r), 'haploid' if args.haploid_genome else 'diploid'] for r in sam_file.references]
    
    if len(regions_list)==0:
        print('\n%s: No valid regions found.'  %str(datetime.datetime.now()), flush=True)
        sys.exit(2)
    
    for contig_data in regions_list:
        if contig_data[0] in ['chrY','Y']:
            contig_data[3]='haploid'
        elif contig_data[0] in ['chrM','M']:
            contig_data[3]='haploid'
        elif contig_data[0] in ['chrX','X']:
            contig_data[3]='haploid' if args.haploid_X else 'diploid'        
    
    tuple_regions_list=[tuple(x) for x in regions_list]
    return tuple_regions_list



def get_chunks(regions_list, cpu, max_chunk_size=500000, min_chunk_size=10000):
    chunks_list=[]
    total_bases=sum(region[2]-region[1]+1 for region in regions_list)
    

    chunksize=min(max_chunk_size, max(min_chunk_size, total_bases//cpu+1))

    for region in regions_list:
        contig=region[0]
        start=region[1]
        end=region[2]
        ploidy=region[3]
        for chunk in range(start, end, chunksize):
            chunks_list.append({'chrom':contig, 'start': chunk, 'end' : min(end, chunk + chunksize), 'ploidy':ploidy})
    
    
    return chunks_list
            
def run_cmd(cmd, verbose=False, output=False, error=False):
    stream=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = stream.communicate()
    
    stdout=stdout.decode('utf-8')
    stderr=stderr.decode('utf-8')
    
    if stderr and error:
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
