#generate_candidate_pileups.py

This module allows us to generate pileup matrices for a given range of reference positions, a BAM and FASTA file. The main function is  generate.

generate(chr_num,start,end,sam_path,fasta_path,out_path)
Inputs:
chr_num - chromosome name (string)
start - starting reference position (int)
end - ending reference position (int)
sam_path - path to BAM file (string)
fasta_path - path to reference FASTA file (string)
out_path - path where the output should be stored (string)

Outputs:
Compressed binary files containing numpy matrices.

This function calls get_candidates() to determine the list of candidate variant sites in the range chr_num:start-end, and stores them in a compressed binary file named 'chr_start_end_candidates.npz'. Then generate() divides the range chr_num:start-end into contiguous 1mb regions. Each region is further divided into ten parts, and pileup information for the parts
are generated in parallel by calling create_pileup() on each 100kb subregion. 

create_pileup(options)
Input:
options - list containing [chr_num,start,end,sam_path,fasta_path], where start and end correspond to the 100kb subregion.

Outputs:
It returns a tuple (pileup, rnames, subregion) where,
pileup - is a numpy matrix of pileups for all the reads in the 100kb subregion, padded by 100 bases on each end to allow us to store a span of bases at the base positions near the ends of 100kb region.
rnames - is an array of read names in the order they appear in pileup matrix
subregion - name identifying the region by its starting and ending point

pileup is matrix of dimension 100200 x N x 5, where N is the number of reads falling in the 100kb(with 100 base padding on each end ) region, including only primary reads.
In the third dimension, first four matrices correspond to presence of a given base at a given reference position, with value 1 if it is a reference match, -1 if non-match, and 0 if that base was not found at that position.

Once the tuple from create_pileup() are returned for each 100kb subregion, they are stored in a compressed binary file named 'chr_start(of the 1mb region)_.npz'.
This file stored each tuple by storing pileup matrix and rnames list for each 100kb region with names 'start-end(of 100kb subregion)' and 'start-end(of 100kb subregion)-rnames' respectively.

The files can be opened and analyzed as follows:
f=np.load(path to .npz file)
f.files() #shows the names of pileup matrix and readname list for each 100kb subregion
f[matrix_name] # returns the matrix with name 'matrix_name' as shown in the command f.files()
