# FFP
Feature Frequency Profile(FFP) two core programs


Prerquirement
A. GCC(g++) version 4.7.1+
B. Google sparse hash library. Can be download here: https://github.com/sparsehash/sparsehash
C. zlib version 1.2.8+. Can be download here: http://www.zlib.net/


1. FFP_compress.cpp
Compile option: g++ -std=c++11 -o (execute name) (this script) -lz
Run example: [Program path][options][input file path][output file path] 

[options(parameters)]
-h  show help, show options
-s  [INT] feature size
-a  take amino acids sequence
-c  convert and accept nucleotide(AGCT) code to RY code
-k  [STR] manual input of word alphabets: For example input 'HJKL' as ['H', 'J', 'K', 'L'] set
-r  disable reverse complement counting
-n  output ratio instead of frequency
-u

-V
-b
-t

[note]
Using [-a], amino acids, automatically turn [-r], disable reverse complement counting, becase peptide sequences have direction(start code -> stop codon), however, nucleotide(genome) sequence actually is double helix that has reverse complement sequence but usually not shown in FASTA file.


[input]
FASTA format peptide or nucleotide sequence files. 


[output]
Data compressed Feature Frequency Profile



2. JSD_matrix.cpp
Compile option: g++ -std=c++11 -pthread -o (execute name) (this script) -lz
Run example: [Program path][options][input files path] > [output file path(standard output)]


[options(parameters)]
-h


[note]

[input]
the output files of 'FFP_compress' which are Feature Frequency Profile(FFP)s

[output]
Standard output of low triangular distance matrix



Limitation:
Generally longer feature length(l-mer) consume more memory.
So far vocabulary size up to 20(amino acids) ^ 24(feature length) was used and available.
