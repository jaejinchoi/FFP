# FFP
Feature Frequency Profile (FFP); two core programs  
A code deposit for "A genome Tree of Life for the Fungi kingdom", JaeJin Choi and Sung-Hou Kim (2017), PNAS.  


## Requirements  
- GCC(g++) version 4.7.1+  
- Google sparse hash library. Look here: https://github.com/sparsehash/sparsehash  
- zlib version 1.2.8+. Look here: http://www.zlib.net/  
  
  
## Tutorial / Supplement
* A tutorial you can walkthrough here: ![Tutorial](example)
* Additional fungi study supplement files (e.g., tree newick and divergence matrix) are here: ![Supplement](fungi_tree_supplement)  
  
  
## Two versions available: Binary, and Text(string)    
* Latest update: 2v.3.0 (2019-2-28); see the ![update history](versions/update_history.txt), list all versions ![versions](versions)  

FFP_binary version; ![FFP_bin, 2v.3.0](versions/2v.3.0/)  
FFP_text version; ![FFP_txt, 2v.2.1](versions/2v.2.1/FFP_txt)  

Old text based FFP 2v.1.0 (before 2018-8); ![2v.1.0](versions/2v.1.0)  
  
## 1. FF Profiler  
Compile: g++ -std=c++11 -o (execute name) (this script) -lz  
Run example: [Program path][options][input file path][output file path]  
When preparing input files save different taxons in separated files  

### [Arguments]
* -h  
    Show options  
* -v  
    Show version 
* -s [INT]  
    Feature size (l-mer)  
* -e [INT]  
    Feature size (l-mer) range end. Functional only when measuring vocabulary size using with -V  
* -a  
    Takes 20 amino acids sequences as inputs  
* -c  
    Convert and accept nucleotide bases (A, G, C, T) into RY bases (Purine, Pyrimidine)
* -k [STR]  
    Manual input of alphabet base string. For example input 'HJKL' is ['H', 'J', 'K', 'L'] bases  
* -r  
    Disable reverse complement accounting. Any bases set other than AGCT code will disable reverse complement accounting  
* -n  
    Output frequency, i.e., feature count / total feature count  
* -u  
    Accept masked (lower confidence; softmasked) letters which are in lower cases in FASTA format  
* -V  
    Count vocabulary size at given range of feature length (l-mer)  
* -b [LONG LONG]  
    Bottom count limit. Remove features, the count below [-b]  
    Default = 1
* -t [LONG LONG]  
    Top count limit. Remove features, the count above [-t]  
    Default = 0 = maximum  
* -B [Double]  
    Bottom entropy limit. Remove features, the string entropy below [-B]  
    Default = 0.0
* -T [Double]  
    Top entropy limit. Remove features, the string entropy above [-T]  
    Default = 1.0 = maximum
    

### [Note]

When using option [-a], amino acids input, [-r] turns on automatically that disable reverse complement accounting, because peptide sequences are single strand that have direction (start code -> stop codon). However, nucleotide sequences input is considered as a double helix, of 'forward' and 'backward' strands, and account reverse compliment as default option.

Use [-r] option that turn off reverse compliment accounting if input is single strand nucleotide sequences such as ribosomal DNAs.

Use [-V] option along with [-s], [-e] and feature filtering arguments to estimate a range of optimal l-mer. In general, use [-b 2], remove any feature count less than 2, to determine a l-mer where vocabulary complexity started to maximize. FF Profiler will determine and stop to a point where vocabular size drops in range of [-s] and [-e], or without giving [-e] (default value is 0) will continue search the point by incrementing the l-mer.  
* Heuristically, an optimal l-mer for amino acids was 13 and for nucleotides was 23 or 24 (but indefinitive). A population's optimal l-mer likely follow the majority optimal l-mer.  

### [Input]
FASTA format peptide or nucleotide sequence files. 


### [Output]
zlib compressed Feature Frequency Profile.



## 2. JSD Calculator  
Compile: g++ -std=c++11 -pthread -o (execute name) (this script) -lz  
Run example: [Program path][options][input files path] > [output file path (standard output)]  

### [Arguments]

* -h  
    Show options  
* -t [INT]  
    Number of threads for multiprocessing  
    Heuristically, an adequate thread number is the number of cpu cores - 1. Default is 5
* -r [PATH]  
    Input previous matrix, and add more items to the matrix without calculating a whole. File name is used as an item name but no longer than 10 characters.
* -d  
    Output Jensen-Shannon distance matrix instead of Jensen-Shannon divergence matrix which is default option.  
    Square root(JS divergence) = JS distance  
* -s  
    Output a symmetric matrix instead of a low triangular matrix which is default output.  

### [Note]
[-r] input low triangular divergence or distance matrix. This requires all pair-wise output of 'FF Profiler'


### [Input]
The output files of 'FFP_compress' which are zlib compressed Feature Frequency Profiles (FFPs)


### [Output]
Standard output of a low triangular Jensen-Shannon divergence or distance matrix in PHYLIP format (not .tsv).  
Support a symmetric matrix output using [-s].

## Limitation
Generally, longer feature lengths (l-mer) consume more memory and time.  
In fungi proteome study the largest proteome has 35,274 proteins containing 10,866,611 amino acids, this program worked for feature length up to 24 amino acids.

