# FFP
Feature Frequency Profile(FFP) two core programs


## Requirements  
A. GCC(g++) version 4.7.1+  
B. Google sparse hash library. Can be download here: https://github.com/sparsehash/sparsehash  
C. zlib version 1.2.8+. Can be download here: http://www.zlib.net/  


## 1. FFP_compress.cpp
Compile option: g++ -std=c++11 -o (execute name) (this script) -lz
Run example: [Program path][options][input file path][output file path] 


### [options(parameters)]

* -h
    Show options  
* -s [INT]  
    Feature size (l-mer)  
* -a
    Takes 20 amino acids sequences as inputs  
* -c
    Convert and accept nucleotide bases into RY bases 
* -k [STR]  
    Manual input of alphabet bases string. For example input 'HJKL' is ['H', 'J', 'K', 'L'] bases  
* -r
    Disable reverse complement counting. Any bases set other than AGCT code will disable reverse complement  
* -n
    Output frequency into ratio. A feature count / Total feature count  
* -u
    Accept masked letters which are lower cases. FASTA format generally expresses masked bases into lower cases  
* -V
    Measure vocabular size at given range of feature length (l-mer)  
* -b [LONG LONG]  
    Bottom limit. Remove any feature count less than -b  
    Default = 1
* -t [LONG LONG]  
    Top limit. Remove any feature count larger than -t  
    Default = 0 = maximum  
    

### [note]
Using [-a], amino acids, automatically turn [-r], disable reverse complement counting, becase peptide sequence have direction(start code -> stop codon), however, nucleotide(genome) sequence actually is double helix that has reverse complement strand.


### [input]
FASTA format peptide or nucleotide sequence files. 


### [output]
Data compressed Feature Frequency Profile




## 2. JSD_matrix.cpp
Compile option: g++ -std=c++11 -pthread -o (execute name) (this script) -lz
Run example: [Program path][options][input files path] > [output file path(standard output)]


### [options(parameters)]

* -h
    Show options  
* -c [INT]
    Specific integer code of single or escape character for delimiter in bewteen a set of [Feature|Value]  
    Default is 13, '\n', linebreak
* -t [INT]
    Number of tread for mulriprocessing  
    Heuristically, an adequate thread number is the number of cpu cores - 1. Default is 5
* -r [PATH]
    Input previous matrix, and add more items to the matrix without calculating a whole    
* -d  
    Output Jensen-Shannon distance matrix instead of Jsensen-Shannon divergence matrix which is default
    
### [note]
[-r] input accept low triangular distance matrix, and it requires all pair-wise FF_Profiles


### [input]
the output files of 'FFP_compress' which are Feature Frequency Profile(FFP)s


### [output]
Standard output of low triangular Jensen-Shannon divergence, or distance, matrix
Square root(JS divergence) = JS distance


## Limitation:
Generally, longer feature lengths(l-mer) consume more memory and time.  
So far vocabulary size up to 20(amino acids) ^ 24(feature length) was used and tested.  
