# generate_candidate_pileups.py

This module allows us to generate pileup matrices for a given range of reference positions, a BAM and FASTA file. The main function is  generate().

## generate(chr_num,start,end,sam_path,fasta_path,out_path)
###### Inputs:
* params - dictionay of paramters
* mode - testing or training mode (string),default=training

the dictionary params has following keys:
1. chr_num - chromosome name (string)
2. start - starting reference position (int)
3. end - ending reference position (int)
4. sam_path - path to BAM file (string)
5. fasta_path - path to reference FASTA file (string)
6. out_path - path where the output should be stored (string)

###### Outputs:
Compressed text file named 'chr_start_end.gz'

generate() function calls get_candidates() to determine the list of candidate variant sites in the range chr_num:start-end. Then, the function generate() divides the range chr_num:start-end into contiguous 1mb regions. Each region is further divided into ten parts, and pileup information for candidate sites in each 100kb subregion is generated in parallel by calling create_training_pileup(). 

Each line of the text file contains all the information pertaining to a candidate site. The format of this text file is follows:
"Position : Genotype : Read 1 name | Read 2 name |... : binary part of pileupmatrix : non-binary part of pileup matrix"


Pileup matrix for a candidate site has dimension=32 x 101 x 8. First dimension corresponds to individual reads. Second dimension correponds to reference base position, with the middle column corresponding to the candidate position. Third dimension corresponds to 8 feature matrices, defined as:
- first four matrices correspond to presence of A, G, T, C base at a given reference position and read, with value 1 if the base is found at that position in the read, and 0 otherwise.
- fifth matrix corresponds to presence of deletion at a given reference position and read, with value 1 if it is a deletion, and 0 otherwise.
- sixth matrix stores reference matches, with value 1 for match and 0 otherwise.
- seventh matrix consistion of all zeros, except ones in the middle column corresponding to the candidate site
- eighth matrix stores base quality information for the base at given reference position and read.

These text files can be read using *read_pileups_from_file(file_path)* function from utils.py

To run this script from commandline use the following syntax:

`python generate_candidate_pileups.py -r chrom:start-end -m training -bam path_to.bam -ref path_to.fa -vcf path_to.vcf -o path_to_output`


