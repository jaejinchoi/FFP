# Example; Tutorial  

## See requirements from main page  
- All proteome data used obtained from NCBI refseqDB, both fungi and non-fungi proteomes as given TAXIDs.  
- Two randoms used as a root made by shuffling the residue order of each peptide sequences in proteome.    

### Arguments  

* Run FF Profiler. Here, using feature length 13 and "FFP_13" as a FF Profiler's output folder  
./FFP_compress -a -s 13 -n 931890 ./FFP_13/931890  
./FFP_compress -a -s 13 -n 332648 ./FFP_13/332648  
./FFP_compress -a -s 13 -n 367775 ./FFP_13/367775  
./FFP_compress -a -s 13 -n 418459 ./FFP_13/418459  
./..  
./..  
./..  
./FFP_compress -a -s 13 -n R990650 ./FFP_13/R990650  


* Run JSD Caculator. Here, using "16_items_13.matrix" as an output  
./JSD_maxtrix -t [number of threads] ./FFP_13/* > 16_items_13.matrix  


* Contstruct a tree using BIONJ and the JS-Divergence matrix



![Workflow](FFP_flowchart3.jpg)
