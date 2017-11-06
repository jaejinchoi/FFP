# Example; Tutorial  

## See requirements from main page  
- 14 proteome data used obtained from NCBI refseqDB which are enlisted in "16_items_list.ods" and "16_items_list.xlsx".  
- 2 randoms used as a root made by shuffling each protein sequence's order in a proteome.  

## Follow  

### 1. Run FFP Profiler  
*-Input is proteome (amino acids), -a
*-Feature length (l-mer) is 13, -s 13
*-Normalize output (frequency), -n
*-Output folder is "./FFP_13"

./FFP_compress -a -s 13 -n 931890 ./FFP_13/931890  
./FFP_compress -a -s 13 -n 332648 ./FFP_13/332648  
./FFP_compress -a -s 13 -n 367775 ./FFP_13/367775  
./FFP_compress -a -s 13 -n 418459 ./FFP_13/418459  
./..  
./..  
./..  
./FFP_compress -a -s 13 -n R990650 ./FFP_13/R990650  


### 2. Run JSD Caculator
*-Using 3 threads, -t 3  
*-Standard output to "16_items_13.matrix", is asymmetric matrix  

./JSD_maxtrix -t 3 ./FFP_13/* > 16_items_13.matrix  

### 3. Construct a tree using Neighbor-Joinging (NJ) method from the divergence matrix 
*-You can use either BIONJ or NJ. However, BIONJ requires to input a symmetric matrix

![Workflow](FFP_flowchart3.jpg)
