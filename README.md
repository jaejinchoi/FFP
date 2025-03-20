# Feature Frequency Profile (FFP)

An alignment-free sequence comparison based on natural language analysis in information theory (e.g., k-mer), which is primary developed to compare whole-genomic sequences such as genomes, proteomes and transcriptome, but it can also compare sequences of English alpabets (books or scripts) and custom alphabet sets other than nucleotide (four letters) or amino acids (20 letters). It constructs reprensenting profiles and determines their distances for further visualization and interpretation.

Please cite one of the publications below if you are utilzing programs provided here for your publication.

* "Whole-proteome tree of life suggests a deep burst of organism diversity", JaeJin Choi and Sung-Hou Kim (2019), PNAS.
* "A genome Tree of Life for the Fungi kingdom", JaeJin Choi and Sung-Hou Kim (2017), PNAS.

## Requirements

Mainly developed and tested in linux environment (Ubuntu 20.04+), and so the command blocks provided below.  
Please contact JaeJin Choi (<jaejinchoi@berkeley.edu>) if you have questions or comments regarding the program codes, bug, or usage.

* GCC (g++) version 4.7.1+  
Any recent g++ versions that support c++11

* Google sparse hash library  
<https://github.com/sparsehash/sparsehash>

```console {ubuntu 24.04}
sudo apt update
sudo apt install libsparsehash-dev
```

* zlib version 1.2.8+  
<https://github.com/madler/zlib> or <http://www.zlib.net/>

```console {ubuntu 24.04}
sudo apt update
sudo apt install zlib1g-dev
```
  
## Tutorial / Supplement

* A tutorial you can walkthrough here: ![Tutorial](example)
* Additional fungi study supplement files (e.g., tree newick and divergence matrix) are here: ![Supplement](fungi_tree_supplement)

## Version compatibility

* The first and the second numberings indicate program compatibility between FF Profile and FFP distance calculation. For instance, any versions between 2v.3.x are compatible but not with any 2v.4.x.
* Incompatibility is due to output file format change during improvement or adding more functions difficult to unify and resolve. Thus, be cautious and when using different versions.  

## Versions

* The latest ![2v.4](versions/2v.4.x)  
* For usage and compiling options, check individual version folder.  
* Old text based ![FFP 2v.1.0 (before 2018-8)](versions/2v.1.0)  
* ![update history](versions/update_history.txt), list all versions ![versions](versions)  

## Note

Typically, longer feature lengths (_ls_) consume more memory. In fungi whole proteome study, the largest proteome has 35,274 proteins containing 10,866,611 amino acids, and the version used worked for feature lengths up to 24 amino acids (l=24).
