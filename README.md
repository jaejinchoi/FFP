# Feature Frequency Profile (FFP)  
Program codes deposit for  
* "Whole-proteome tree of life suggests a deep burst of organism diversity", JaeJin Choi and Sung-Hou Kim (2019), PNAS.  
* "A genome Tree of Life for the Fungi kingdom", JaeJin Choi and Sung-Hou Kim (2017), PNAS.  

Please cite one of the publications above.  

## Requirements  
- GCC(g++) version 4.7.1+  
- Google sparse hash library. Look here: https://github.com/sparsehash/sparsehash  
- zlib version 1.2.8+. Look here: http://www.zlib.net/  

Welcome to have questions or opinions regarding the program codes, bug, and usage. 
* Please contact JaeJin Choi (jaejinchoi@berkeley.edu) 
  
## Tutorial / Supplement
* A tutorial you can walkthrough here: ![Tutorial](example)
* Additional fungi study supplement files (e.g., tree newick and divergence matrix) are here: ![Supplement](fungi_tree_supplement)  
  
## Version compatibility
* The first and the second numberings indicate program compatibility between FF Profile and FFP distance calculate. For instance, any versions between 2v.3.x are compatible but not with any 2v.4.x.  
* Incompatibility is due to output file format change during improvement or adding more functions which is difficult to resolve. Thus, please, beware using different versions.  

## Versions
* The latest ![2v.4](versions/2v.4.x)  
* For usage and compiling options, check individual version folder.  
* Old text based ![FFP 2v.1.0 (before 2018-8)](versions/2v.1.0)  
* ![update history](versions/update_history.txt), list all versions ![versions](versions)  


### [Note]
Typically, longer feature lengths (ls) consume more memory. In fungi whole proteome study, the largest proteome has 35,274 proteins containing 10,866,611 amino acids, and the version used worked for feature lengths up to 24 amino acids (l=24).  

