# Feature Frequency Profile (FFP) wrapper script

Automate iteration of feature lengths and determine optimum feature length (_l_) range.

## <span style="color: green;">Requirements</span>

This Python script is largely dependent on standard packages need not to install; however, need some extra packages for core functions. Required packages can be easily find and install using PIP.

* Python 3.3+ _# check the version with --version_
* Pandas _# summarize results_
* Matplotlib _# plotting using pandas.dataframe()_
* Dendropy _# build and compare trees_

```console
pip3 install pandas metplotlib dendropy
```

## <span style="color: green;">Copy of help page</span>

Usage: [program_path] [options] [path_to_input_files]

[options]

* -h or --help, show help
* -v or --version, show version and citation information

### Feature length iteration from [-s] to [-e]

* -s [int], from feature length (default=1)
* -e [int], to feature length (default=1)

### Process behavior

* -S [working_folder_path], a whole process and all output are constriucted in the designated folder
* -w or --overwrite, overwrite all output (will remove old saved files and folders)
* -f, force to search up to given feature lengths, [-e], beyond upperboud  
Upperbound is a point where the first maximum distance (between ingroup taxons) found
* -t or --thread [int], a number of threads or cores to use
* --outgroup [str], designate outgroup taxon; multiple outgroups separated by comma. Specified outgroup(s) will be used for rooting trees and determining upperbound (of ingroups)

### (external) Program paths and arguments

#### 1. Construct FF-Profile

* --ffp_ppath [path], the path to FF-Profiling program
* --ffp_parg [str], arguments for FF-Profiling program except input and output

#### 2. Construct (pairwise) Distance Matrix

* --dmat_ppath [path], the path to dist_matrix program
* --dmat_parg [str], arguments for dist_matrix program except input and output

#### 3. Construct Tree

* --tree_ppath [path], the path to (external) tree build program of choice
Without --tree_ppath will use dendropy for tree construction using --tree_method
* --tree_method [str], tree construction method (default=nj); currently support ['nj', 'upgma']
* --tree_parg [str], arguments for (external) tree build program except input and output
Will call external tree build program: [tree_ppath] [tree_parg] [input] [output], in this order

* --root_method [str], rooting by "node" or "edge" (default=node)
* --tree_schema [str], tree format types: Newick, NEXUS, NeXML (default=Newick)

#### 4. Summary and Plot

--rf_norm, normalize RF distance; recommended to use (True); default is False
