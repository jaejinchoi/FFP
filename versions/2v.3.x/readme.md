### 2v.3.x  

## [Change]

* From version 2v.3.0, only binary FFP will be continued and supported (for text-based FFP the latest version is 2v.2.1)
* FFP_bin's value presentation is Long Long (for feature count) or Double (for feature frequency).  
  A variable byte definition dependent on platforms because it is byte based.
* JSD_matrix_bin's value presentation is up to 8 significant figures below a decimal point (%.8g)  
* Option [-v] shows a program profile  
* 2v.3.0 fully supports user-defined option [-k] 

2022-11-5

* An option '-p' to print features and the count in plain format, which takes more space

2022-4-5

* FFP distance; add option to use TAB a row name separater, instead of PHYLIP format that limit row names up to 9 characters

2020-2-24

* FFP distance; a valid-time of 'q_f_buf' may cause adverse consequences and so there was a code rearrangement.  

2020-1-13

* Unadvised to input large genome files. It shown to loss/malfunction during zlib compression.  

## [Limitation]

By the design the maximum feature length (l-mer) is supported up to  
1785 for Purine-Pyrimidine (RY) encoded sequences (2 letters),  
892 for Nucleotide sequences (4 letters) and 0,1,2 genotype sequences (3 letters), and  
357 for Amino acids sequences (20 letters).  

The maximum l-mer is calculated by two steps, include a number of customized letters (parameter -k).  
a. Estimate bits_per_letter, for instance, 2 bits (2^2=4) required for 4 letters and 5 bits (2^5=32) required for 20 letters  
b. floor(1,785 / bits_per_letter) >= your maximum l-mer  


## FF Profile(FFP); FFP_x.cpp

**Compile:**

```console
g++ -std=c++11 -o FFP_bin_2v.3 FFP_bin_2v.3.cpp -lz
```

Run example: [Program path][arguments][input file path][output file path]  
Each input file represent one operational taxon unit (OTU) in a tree.  

### [Input]

FASTA format sequence files.  

### [Output]

zlib compressed Feature Frequency Profile.  

### [Arguments]

* -h  
    Show help  
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
    Default is off; e.g., AAA -> 122, then the reverse complement is TTT -> 122 but not shown  
* -p  
    Print plain features and the counts without binary compression  
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

Use [-r] option to turn off reverse compliment accounting if input is single strand nucleotide sequences of direction such as ribosomal RNAs.

Use [-V] option along with [-s], [-e] and feature filtering arguments to estimate a range of optimal l-mer. In general, use [-b 2], remove any feature count less than 2, to determine a l-mer where vocabulary complexity started to maximize. FF Profiler will determine and stop to a point where vocabular size drops in range of [-s] and [-e], or without giving [-e] (default value is 0) will continue search the point by incrementing the l-mer.  

## FFP distance calculate; JSD_matrix_x.cpp

**Compile:**

```console
g++ -std=c++11 -pthread -o JSD_matrix_bin.2v.3 JSD_matrix_bin.2v.cpp -lz
```

Run example: [Program path][arguments][input files path] > [output file path (standard output)]  

### [Input]

zlib compressed Feature Frequency Profiles (FFPs).  

### [Output]

Standard output of a low triangular distance matrix.  

### [Arguments]

* -h  
    Show options
* -v  
    Show version
* -t [INT]  
    Number of threads for multiprocessing. An adequate thread number is a total cpu thread - 1. Default is set to 5.    
* -r [PATH]  
    Input previous matrix and append more items. Should match previous matrix format (e.g., PHYLIP or tab)  
* -d  
    Convert Jensen-Shannon Divergence to Jensen-Shannon Distance, which is equivalent to square root of Jesen-Shannon Divergence.
* -T 
    Use TAB as a separator between row names and distances. Default is PHYLIP format that limit row names up to 9 characters.  
* -s  
    Output a symmetric matrix instead of a low triangular matrix which is default output.  

### [Note]

[-r] input previously generated low triangular divergence or distance matrix. This requires all pair-wise output of 'FF Profiler'. Aware to specify "TAB" delimited and "PHYLIP" formatted distance matrices.  
