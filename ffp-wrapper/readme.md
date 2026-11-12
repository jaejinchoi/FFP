# Feature Frequency Profile (FFP) wrapper script

Automate iteration of feature lengths, from FF-Profiling to calculating distance matrices, and then summarize Robinson-Foulds (RF) distances and the upperbound to help determine the optimum feature length (_l_) range.

The process is divided into three subcommands, connected through formatted output paths:

1. **profile**: construct FF-Profiles for each _l_
2. **distance**: calculate (pairwise) distance matrices for each _l_
3. **optimum**: build trees from the distance matrices, calculate RF distances between trees of adjacent _l_s, and then create a summary and plot

## <span style="color: green;">Requirements</span>

This Python script is largely dependent on standard packages need not to install; however, need some extra packages for core functions. Required packages can be easily find and install using PIP.

* Python 3.12+ _# check the version with --version_
* Pandas _# read distance matrices and summarize results_
* Matplotlib _# plotting using pandas.dataframe()_
* Dendropy _# build and compare trees_

```console
pip3 install pandas matplotlib dendropy
```

Compiled FFP programs are also required: the FF-Profiling program (e.g., FFP_bin_2v.4) and the distance matrix program (e.g., JSD_matrix_bin_2v.4.1).

## <span style="color: green;">Workflow and output</span>

Run **profile**, then **distance** with the FF-Profile folders, then **optimum** with the distance matrices. **distance** and **optimum** accept inputs given by users as long as the names match the format, [prefix].[l] (e.g., fp.10 and matrix.jsdiv.10); inputs are sorted by _l_ (the numeric name suffix) in ascending order.

```console
python3 ffp_wrapper.py profile -S [working_folder] -s 3 -e 20 --ppath [path_to]/FFP_bin_2v.4 --parg "-n" [path_to_sequence_files]

python3 ffp_wrapper.py distance -S [working_folder] --ppath [path_to]/JSD_matrix_bin_2v.4.1 --parg "-d jsdiv" --outgroup [outgroup_taxon] [working_folder]/fp.*

python3 ffp_wrapper.py optimum -S [working_folder] --outgroup [outgroup_taxon] --rf_norm [working_folder]/distance_matrix/matrix.jsdiv.*
```

Output structure in the working folder [-S]:

```console
[working_folder]/
├── fp.[l]/                                        # 1. profile: FF-Profiles of input sequences for each l
├── distance_matrix/matrix.[distance type].[l]     # 2. distance: e.g., matrix.jsdiv.10
├── tree/tree.[distance type].[l]                  # 3. optimum: (rooted) trees, e.g., tree.jsdiv.10
├── summary.csv                                    # 4. optimum: upperbound, minimum RF distance, and per l summary
├── search_optimum.png                             # 4. optimum: RF distance plot with the upperbound (red dashed line)
├── [date]_out.txt                                 # standard output of the called programs
└── [date]_error.txt                               # standard error of the called programs
```

* Upperbound is a point (_l_) where the first maximum distance between ingroup taxons (excluding outgroups) found
* RF distance at _l_ is between trees of _l_ and the next _l_ (in ascending order) among given distance matrices
* The optimum _l_ range is to be determined by a user, using the RF distance plot and upperbound; the minimum RF distance in summary.csv is a guide

## <span style="color: green;">Copy of help page</span>

Usage: python3 ffp_wrapper.py [-h] [-v] {profile,distance,optimum} [options]

* -h or --help, show help; [subcommand] -h shows help of each subcommand
* -v or --version, show version and citation information

### Options shared by subcommands

#### Feature length iteration from [-s] to [-e] (profile)

* -s [int], from feature length (default=1)
* -e [int], to feature length (default=1)

#### Process behavior (all subcommands)

* -S [working_folder_path], a whole process and all output are constructed in the designated folder
* -w or --overwrite, overwrite existing output files of the same path; without it, existing (non-empty) output files are skipped. Folders and other files are kept
* -t or --thread [int], a number of threads or cores to use (default=available cores - 1)

### 1. profile

Construct FF-Profiles of input files for each feature length; saved to [-S]/fp.[l]

Usage: python3 ffp_wrapper.py profile [options] --ppath [path] [path_to_input_files]

* --ppath [path], the path to FF-Profiling program (required)
* --parg [str], arguments for FF-Profiling program except feature length, input and output  
Feature length is given as -s [l] for each _l_ from [-s] to [-e]; -s in --parg overrides it (warned)
* [path_to_input_files], input sequence files

### 2. distance

Construct (pairwise) distance matrix from FF-Profile folders made by profile, sorted by feature length (the numeric folder name suffix, [l]) in ascending order; saved to [-S]/distance_matrix/matrix.[distance type].[l]

Usage: python3 ffp_wrapper.py distance [options] --ppath [path] [path_to_profile_folders]

* [path_to_profile_folders], FF-Profile folders, e.g., [-S]/fp.*

* --ppath [path], the path to dist_matrix program (required)
* --parg [str], arguments for dist_matrix program except input and output  
-t [thread] and -T (tab delimited) are added if missing. Distance type is taken from -d in --parg (default=jsdiv)
* -f or --force, force to search all given feature lengths beyond upperbound  
Without -f, stop at the upperbound
* --outgroup [str], designate outgroup taxon; multiple outgroups separated by comma  
Specified outgroup(s) are excluded when determining upperbound (of ingroups)

### 3. optimum

Construct (rooted) trees from distance matrices made by distance, sorted by feature length (the numeric file name suffix, [l]) in ascending order; calculate RF distance between trees of adjacent _l_s, then summary and plot. Trees are saved to [-S]/tree/tree.[distance type].[l]

Usage: python3 ffp_wrapper.py optimum [options] [path_to_distance_matrices]

* [path_to_distance_matrices], distance matrix files, e.g., [-S]/distance_matrix/matrix.jsdiv.*  
Lower triangular or square matrix
* --outgroup [str], designate outgroup taxon; multiple outgroups separated by comma  
Specified outgroup(s) will be used for rooting trees and excluded when determining upperbound (of ingroups); included in tree building and RF distance

#### Construct tree

* --method [str], dendropy tree construction method (default=nj); currently support nj and upgma
* --root_method [str], rooting by "node" or "edge" (default=node)
* --schema [str], tree format types: newick, nexus, nexml (default=newick)
* --ppath [path], the path to (external) tree build program of choice  
Without --ppath will use dendropy for tree construction using --method
* --parg [str], arguments for (external) tree build program except input and output  
Will call external tree build program: [ppath] [parg] [input] [output], in this order

#### Summary and Plot

* --rf_norm, normalize RF distance; recommended to use (True); default is False
