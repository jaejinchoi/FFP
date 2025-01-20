# python3.3 <=


### README
# You can utilize this python wrapper to conveniently iterate over feature lengths and to optimize through resulting plots.
# Also, you can modify/distribute, at your risk, to run this script for your interest
###

import os, sys, shutil, getopt, io
import subprocess, shlex, time

from time import gmtime, strftime

from multiprocessing import cpu_count, Pool
from typing import overload

import dendropy
from dendropy.calculate import treecompare

import pandas as pd
import numpy as np
import statistics

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
#from Bio.Phylo.TreeConstruction import DistanceCalculator


#DEFAULT SET
DEFAULT_SAVE_FFP_FOLDER_PREFIX=""
DEFAULT_SAVE_MATRIX_FOLDER_NAME=""
DEFAULT_SAVE_TREE_FOLDER_NAME=""
DEFAULT_DEDROPY_TREE_METHOD_LIST=["nj", "upgma"]

#2012.9.28, act as qsub-cluster function(limit core number- input job) using multiprocessing
#2017-5 using Pool
class aux_core:

    def create_folder(
        self
        , folder_path=""
        , overwrite_flag=False
        ):

        if (os.path.exists(folder_path)==True and overwrite_flag==True):
            shutil.rmtree(folder_path)

        if (os.path.exists(folder_path)==False):
            os.mkdir(folder_path)

        print(f"Create working folder: {folder_path} (overwrite={overwrite_flag})")


    def run_multiprocess(
        self
        , proc_name
        , args_list=[]
        , n_cpu=0
        ):

        proc_pool = Pool(processes=n_cpu) #define pool

        ret_list = proc_pool.starmap(proc_name, args_list)

        proc_pool.close()
        proc_pool.join()

        for index, n_item in enumerate(ret_list):
            if (n_item=="FALSE"): #fail to retreive or process result
                print(f"ERROR returning {args_list[index]}")

        if ("FALSE" in ret_list):
            raise ValueError("#HALT: unable to proceed next steps")

        #return ret_list


    def subprocess_call(
        self
        , work_term=''
        , job_output_path=''
        , job_error_path=''
        , cwd_env={}
        ):

        out_file = open(job_output_path, 'w')
        err_file = open(job_error_path, 'a')
        subprocess.call(shlex.split(work_term), stdout = out_file, stderr = err_file, env=cwd_env, shell=False) #give a file handle; shell=False takes arguments as a list
        err_file.close()
        out_file.close()

        return(0)


class ffp_core: #only generate arguments and do the rest out of the function
    def profiling(
        self
        , program_path=""
        , ext_args_str=""
        , load_path_list=[]
        , save_folder_path=""
        , overwrite_flag=False
        , job_output_path=''
        , job_error_path=''
        ):

        args_list=[]

        for load_path in load_path_list:
            save_path = "%s/%s" % (save_folder_path, os.path.basename(load_path))
            work_term = "%s %s %s %s" % (program_path, ext_args_str, load_path, save_path)

            if (overwrite_flag==True or (os.path.exists(load_path) and os.path.getsize(load_path)!=0) ):    
                args_list.append([work_term, job_output_path, job_error_path, cwd_env])
            else:
                print("SKIP: %s" % (work_term))


    def pairwise_dmatrix(
        self
        , program_path=""
        , ext_args_str=""
        , load_folder_path=""
        , save_path='' #, save_folder_path=""
        , cwd_env={}
        , overwrite_flag=False
        , job_output_path=''
        , job_error_path=''
        ):

        file_list=[]

        for root, dirs, files in os.walk(load_folder_path): #recursive search for (ffprofile) files
            for name in files:
                file_list.append("%s/%s" % (root, name))

        work_term = "%s %s %s > %s" % (program_path, ext_args_str, ", ".join(file_list), save_path)

        if (overwrite_flag==True or (os.path.exists(load_path) and os.path.getsize(load_path)!=0) ):    
            args_list.append([work_term, job_output_path, job_error_path, cwd_env])
        else:
            print("SKIP: %s" % (work_term))

        args_list.append([work_term, job_output_path, job_error_path, cwd_env])



class tree_handle:
    # https://dendropy.readthedocs.io/en/main/library/phylogeneticdistance.html
    # https://dendropy.readthedocs.io/en/main/primer/phylogenetic_distances.html
    def dendropy_DistMatrix_to_Tree(
        self
        , load_path=""
        , save_path=""
        , tree_method=""
        ):
        #output: Bool

        taxon_cnt=0
        tree = ""

        pdm = dendropy.PhylogeneticDistanceMatrix()
        str_io = io.StringIO()

        with open(load_path) as read_f: #cannot accept low triangular matrix
            taxon_cnt = int(read_f.readline())
            #print(taxon_cnt)
            str_io.write(read_f.read().strip()) #to remove any blank lines at the end; else, will raise error as proceed
            str_io.seek(0)

            # require to remove any trailing blank lines at the end
            pdm = dendropy.PhylogeneticDistanceMatrix.from_csv( #considering only the upper righ triangle portion
                src=str_io
                , is_first_row_column_names=False
                , is_first_column_row_names=True #taxon names enlisted on the first column
                , is_allow_new_taxa=True #make new
                , delimiter="\t" #reason distance matrix should be tab or comma delimited (unsupport delimit_whitespace)
                #, delimit_whitespaces=True
                )

        str_io.truncate(0)

        if (tree_method=="nj"):
            tree = pdm.nj_tree()
            #print(tree.as_string("newick"))

        elif (tree_method=="upgma"):
            tree = pdm.upgma_tree()
            #print(tree.as_string("newick"))

        else:
            print("Unsupport: %s" % (tree_method))
            return(False)

        if (save_path!=""):
            with open(save_path, "w") as write_f:
                write_f.write(tree.as_string("newick"))


    def external_build_Tree(
        self
        , load_path=""
        , save_path=""
        , argument_str=""
        , overwrite_flag=False
        , cwd_env={}
        ):

        if (os.path.exists(save_path)==True and overwrite_flag==False):
            pass

        else:    
            aux_core().subprocess_call(
                work_term = f"{argument_str} {load_path} {save_path}".strip()
                , job_output_path = job_output_path
                , job_error_path = job_error_path
                , cwd_env = cwd_env
            )


    def reroot_tree( #simplify the code
        self
        , tree_load_path=''
        , tree_save_path=''
        , mrca_clade_list=[]
        , rooting_method=""
        ): #using dendropy

        tree = dendropy.Tree.get_from_path(
            tree_path = tree_load_path
            , schema = "newick"
        )

        tree.encode_bipartitions()
        #tree.is_rooted = True

        #https://dendropy.readthedocs.io/en/v3.12.1/tutorial/treemanips.html
        #reroot options: edge, node, and outgroup (single)
        if (len(mrca_clade_list)>1):

            mrca = tree.mrca(taxon_labels=mrca_clade_list) #mrca = most recent common ancestor of multiple outgroups
            
            if (rooting_method=="node"):
                tree.reroot_at_node(mrca, update_bipartitions=True)
            
            elif (rooting_method=="edge"):
                half_edge = mrca.edge_length / 2.0
                tree.reroot_at_edge(mrca.edge, length1=half_edge, length2=half_edge, update_bipartitions=True, suppress_unifurcations=False)

            else:
                raise ValueError(f"unsupported rooting method: {rooting_method} using {mrca_clade_list}")
            ## midpoint rooting unsupported for multiple outgroups

        else: #pull one outgroup, then halved its edge
            outgroup_node = tree.find_node_with_taxon_label(reroot_clade_list[0])

            if (rooting_method=="node"):
                tree.to_outgroup_position(outgroup_node, update_bipartitions=True)
            
            elif (rooting_method=="edge"):
                half_edge = outgroup_node.edge_length / 2.0
                tree.reroot_at_edge(outgroup_node.edge, length1=half_edge, length2=half_edge, update_bipartitions=True, suppress_unifurcations=False)

            else:
                raise ValueError(f"unsupported rooting method: {rooting_method} using {mrca_clade_list}")

        tree.is_rooted = True

        with open(tree_save_path, "w") as write_f:
            write_f.write(tree.as_string("newick"))


    ## short coming of Robinson-Foulds distance
    # https://ms609.github.io/TreeDist/articles/Generalized-RF.html
    def RF_distance(
        self
        , tree_path_a=''
        , tree_path_b=''
        , tree_format='newick'
        , normalize_flag=False
        ):

        tree_ns = dendropy.TaxonNamespace()
        tree_a = dendropy.Tree.get_from_path(tree_path_a, tree_format, taxon_namespace=tree_ns)
        tree_b = dendropy.Tree.get_from_path(tree_path_b, tree_format, taxon_namespace=tree_ns)

        tree_a.encode_bipartitions()
        tree_b.encode_bipartitions()

        ###tree compare manual
        #https://dendropy.org/library/treecompare.html
        rf_distance = treecompare.symmetric_difference(tree_a, tree_b, False) #unweighed

        if (normalize_flag==True):
            #https://www.rdocumentation.org/packages/phangorn/versions/2.5.5/topics/treedist
            #normalized robinson-fould distance; 2n - 6 is the maximum difference if unrooted and binary trees. n is a number of tips or OTUs
            return float(rf_distance) / float(len(tree_a.taxon_namespace) + len(tree_b.taxon_namespace) - 6)

        else:
            return rf_distance



class do_plot:
    def search_bound(
        self
        , load_path_list=[]
        ):

        pass

    def rb_distance(
        self
        , load_path_list=[]
        ):

        pass


def root_determine(pd_df = pd.DataFrame()):
    ## potential outgroup for rooting purpose. Select the OTU(s) far from the rest.
    rsum = pd_df.sum(axis=1) #so far result the same using mean - std
    root_list = list(pd_df.index[rsum==max(rsum)])
    rsum_mean = statistics.mean(rsum)

    print("\toutgroup/root determine")
    print("Maximum row sum = %.2f, mean row sum: %.2f. suggested root(s) %s" % (max(rsum), rsum_mean, ",".join(root_list)))
    
    return(root_list)






def show_help():
    print('Parameter [option][load_path]')
    print('[load_path] can be either files or directories')
    print("\t#[Set environment]")

    print('-S [folder_path], non-default working folder where all intermediate files and results are saved')

    print('--ffp [path], the path to FF-Profiling program')
    print('--ffp_xarg [path], arguments for ffp_program except input and output')
    
    print('--mat [path], the path to dist_matrix program') 
    print('--mat_xarg [str], arguments for dist_matrix program except input and output')


    print("\t#[feature length arguments]")
    print('-s [int], start feature length(default=1)')
    print('-e [int], end feature legnth(default=0)')
    print('\t-s 0 -e e, without lmer-option, single run')
    print('\t-s > 0 -e 0, single cycle of -s')
    print('\t-s > 0 -e > -s, range work  from -s to -e')


    print('\n#[Set optimum search arguments]')
    print('--rf_dev [float], topological variation deviation upper bound (default=20%?)')
    print('--metric_ubound [int], metric upperbound, stop when any distances reach the metric maximum (information unrelated)')
    print('\tReaching a metric maximum means two information are unrelated, thus, violate any pairs of OTUs are evolutionary unrelated (no common ancestor)')

    print('\n#[Accessory arguments]')
    print('-T, test print for work_term')
    print('-m [int], a number of thread(default = 1)')
    print('\t-1, use maximum cores - 1, when to maximize core use')
    sys.exit()



def show_version():
    #print '2013.8.26 job_submit script designed for universial purpose'
    print('2017.5.19 add more options, and simplifed functions')
    #print 'base on, mk_wrapper2_reserve.py'
    print('support workstation(multiprocessing) and cluster(using SGE, qsub)')
    print('2018.6.24, add qsub (cluster) core usage limit option')

    sys.exit()


if __name__=='__main__':

    try:
        opts, args = getopt.getopt(
            sys.argv[1:]
            , 'ht:s:e:S:m:rR:M:H:Q:'
            , [
                "help"
                , "threads=" #how many cores or threads to use
                , "overwrite" #overwrite all output (will remove old ones); default=False
                # for construct FFP
                , "ffp_path="
                , "ffp_xarg="
                # for consturct distance matrix
                , "dmat_path="
                , "dmat_xarg="

                , "rf_dev="
                , "metric_ubound="

                , "tree_ext="
                ,"tree_xarg="
                , "tree_method="
                , "rf_norm"
                , "outgroup="
                , "root_method="
                ]
            )

    except:
        show_help()

    #ext_arguments_term='' #for external program
    #initialize arguments or flags
    args_dict={}
    args_dict["-w"] = False #do not overwrite
    args_dict["-T"] = False #test run (to check arguments input)
    args_dict["--tree_method"]="nj" #default=nj, neighbor-joining
    args_dict["--rf_norm"]=False #normalize RF by the maximum RF. Maximum RF = 2 * (a number of OTUs) - 6

    #necessary parameters
    l_begin=1 #0 as minimum
    l_end=0 #default=0

    n_cores = 1 #default = 1
    rf_diff_lthreshold = 0.2

    ## default working paths
    save_work_path=os.getcwd()
    ffp_program_path=""
    dist_program_path=""
    tree_program_path=""
    
    #tree manipulation
    reroot_clade_list=[] #one or multiple OTUs for rooting purpose
    root_method="node"

    overwrite_flag=False

    for opt, arg in opts:

        if (opt=="-h" or "--help"):
            show_help()

        elif (opt=="-v"):
            show_version()

        elif (opt=="-t" or opt=="--threads"):
            n_cores = max(1, int(arg)) # least 1

        elif (opt=="-s"):
            l_begin=int(arg)

        elif (opt=="-e"):
            l_end=int(arg)

        elif (opt=="-f" or opt=="--force"):
            l_end=int(arg)

        elif (opt=="--overwrite"):
            overwrite_flag=True

        elif (opt=="--ffp"): #ffprofile program path
            ffp_program_path = os.path.abspath(arg)
            #args_dict["--ffp"] = os.path.abspath(arg)

        elif (opt=="--ffp_xarg"): #arguments for ffprofile program
            args_dict["--ffp_xarg"] = str(arg).strip() #.strip() #external arguments

        elif (opt=="--mat"): #dist matrix program path
            dist_program_path = os.path.abspath(arg)
            #args_dict["--mat"] = os.path.abspath(arg)

        elif (opt=="--mat_xarg"): #arguments for ffprofile program
            args_dict["--mat_xarg"] = str(arg).strip() #.strip() #external arguments

        elif (opt=="--tree_ext"): #dist matrix program path
            tree_program_path = os.path.abspath(arg)
            #args_dict["--tree"] = os.path.abspath(arg)

        elif (opt=="--tree_xarg"): #arguments for ffprofile program
            args_dict["--tree_xarg"] = str(arg).strip() #.strip() #external arguments

        elif (opt=="--tree_method"): #arguments for ffprofile program
            args_dict["--tree_method"] = str(arg).strip().lower() #.strip() #external arguments

        elif (opt=="-S"): #non-default working directory
            save_work_path=os.path.abspath(arg)

        elif (opt=="-m"):
            n_cpu = int(arg)

            if (n_cpu==-1):
                n_cpu = cpu_count()-1 #use this not knowing a core count

        elif (opt=="-T"):
            test_print_flag=True

        elif (opt=="-w"): #overwrite all
            #args_dict["-w"]=True
            overwrite_flag=True

        elif (opt=="--res_mat"): #input previously made distance matrix (add new items)
            args_dict["--res_mat"] = os.path.abspath(arg)

        elif (opt=="--outgroup"): #designate outgroup for rooting purpose
            args_dict["--outgroup"] = str(arg).split(",")

        elif (opt=="--root_methodp"): #designate outgroup for rooting purpose
            args_dict["--root_method"] = str(arg)

    #set derived arguments
    args_dict["prefix_fp_path"] = "%s/fp." % (save_work_path)
    args_dict["prefix_mat_path"] = "%s/distance_matrix/mat." % (save_work_path)
    args_dict["prefix_tree_path"] = "%s/tree/tree." % (save_work_path)


    #pre-check
    if (l_begin>0 and l_end==0): #single l run
        l_end = l_begin

    elif (l_begin<0):
        print("l_start and l_end should be larger than 0, in case of range work")
        sys.exit()


    load_path_list=[]
    instruction_term_list=[]

    argument_term=''

    #qsub_out
    current_time_str = strftime("%Y.%m.%d", gmtime())
    job_output_path="%s/%s_out" % (os.getcwd(), current_time_str)
    job_error_path="%s/%s_error" % (os.getcwd(), current_time_str) #can be avoid by using qsub -j y

    ffprofile_program_path=""
    distmatrix_program_path=""
    work_dir_path=os.getcwd()

    dist_mat_program_path=""
    tree_program_path=""

    rf_distance_list=[]

    cwd_env = os.environ #currento working environment for shell execution

    if (len(args)>0):

        load_path_list=[] # of FASTA formatted sequence files

        for n_path in args: #can be either files or folders
            load_path = os.path.abspath(n_path)
            #prevent file path interruption in linux environment
            #load_path = load_path.replace('(', '\(')
            #load_path = load_path.replace(')', '\)')
            #load_path = load_path.replace(' ', '\ ')
            load_path_list.append(load_path)

        ## 0. Check available programs
        if (os.path.exists(ffprofile_program_path)==False
            or os.path.exists(distmatrix_program_path)==False
            ):
            print(f"Either cores programs {ffprofile_program_path} or {distmatrix_program_path} not exists")
            raise ValueError(f"# Please check if they are compiled and executable")

        # The main iteration per feature length
        for feature_length in range(l_begin, l_end+1):
            output_path_ffp = f"{work_dir_path}/fp.{feature_length}"
            output_path_dmatrix = f"{work_dir_path}/matrix/matrix.fp.{feature_length}"
            output_path_tree = f"{work_dir_path}/tree/tree.fp.{feature_length}"

            ## 1. Construct Feature Frequency Profiles, takes the argument --ffp_xarg
            aux_core().create_folder(
                folder_path = output_path_ffp
                , overwrite_flag=overwrite_flag
            )

            args_list = ffp_core().profiling(
                program_path=ffp_program_path
                , ext_args_str=args_dict["--ffp_xarg"]
                , load_path_list=load_path_list
                , save_folder_path=output_path_ffp
                , n_cores=n_cores
                , cwd_env=cwd_env
                , overwrite_flag=overwrite_flag
                , job_output_path=job_output_path
                , job_error_path=job_error_path
                )

            aux_core().run_multiprocess(
                aux_core().subprocess_call
                , args_list=args_list
                , n_cpu=n_cores
                )

            ## 2. Calculate pairwise distance matrix
            aux_core().create_folder(
                folder_path = output_path_dmatrix
                , overwrite_flag=overwrite_flag
            )

            args_list = ffp_core().pairwise_dmatrix(
                program_path=dist_mat_program_path
                , ext_args_str=args_dict["--mat_xarg"]
                , load_folder_path=output_path_ffp
                , save_path=output_path_dmatrix
                , cwd_env=cwd_env
                , overwrite_flag=overwrite_flag
                , job_output_path=job_output_path
                , job_error_path=job_error_path
                )

            aux_core().run_multiprocess(
                aux_core().subprocess_call
                , args_list=args_list
                , n_cpu=1 #multithreaded in the called program
                )


            ## 3. Build a tree from a distance matrix, takes the argument --dmatrix_xarg
            if (args_dict["--tree_method"] in DEFAULT_DEDROPY_TREE_METHOD_LIST and tree_program_path==""):
                tree_handle().dendropy_DistMatrix_to_Tree(
                    load_path=output_path_dmatrix
                    , save_path=output_path_tree
                    , tree_method=args_dict["--tree_method"]
                    )
            
            ### 3.1 use external tree building program
            elif (tree_program_path!="" and os.path.exists(tree_program_path)):
                tree_handle().external_build_Tree(
                    load_path=output_path_dmatrix
                    , save_path=output_path_tree
                    , argument_str=f"{tree_program_path} {args_dict['--tree_xarg'].strip()} "
                    , overwrite_flag=overwrite_flag
                    , cwd_env=cwd_env
                    )

            else:
                print(f"# Dendropy tree method: {args_dict['--tree_method']} not supported")
                print(f"## Currently supports {DEFAULT_DEDROPY_TREE_METHOD_LIST}")
                print(f"## Or check the external tree program path: {tree_program_path}")
                sys.exit()

            ## 3.1 Reroot the tree if instructed, takes the argument --root
            rerot

        # Create plots for optimization, using pandas

        ## 1. Search bound: 
        ### A. Lowerbound: a point of (matrix) correlation coefficient begin to converge
        ### B. Upperbound: a point of the first maximum distance appears

        ## 2. Robinson-Foulds metric based tree distance trend


        else:
            print("check FF-Profiling program path: %s" % (ffp_program_path))
            print("or check dist matrix program path: %s" % (dist_program_path))

    else:
        print("without input files")

