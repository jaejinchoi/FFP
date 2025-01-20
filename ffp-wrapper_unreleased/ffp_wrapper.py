#python2 and python3 compatible
import os, sys, shutil, getopt, io
import subprocess, shlex, time

from time import gmtime, strftime

from multiprocessing import cpu_count, Pool
from typing import overload

import dendropy
from dendropy.calculate import treecompare

import pandas as pd
import statistics

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
#from Bio.Phylo.TreeConstruction import DistanceCalculator

#2012.9.28, act as qsub-cluster function(limit core number- input job) using multiprocessing
#2017-5 using Pool
def spontaneous_process(proc_name, args_list=[], n_cpu=0):

    proc_pool = Pool(processes=n_cpu) #define pool

    ret_list = [proc_pool.apply_async(proc_name, args=n_arg) for n_arg in args_list]

    proc_pool.close()
    proc_pool.join()

    result_list=[]

    for n_item in ret_list:
        result_list.append(n_item.get()) #[] = failed, []!= success


    del proc_pool #refresh(remove) after use

    return result_list


        
def multiprocess_unit(work_term='', job_output_path='', job_error_path='', cwd_env={}):

    out_file = open(job_output_path, 'w')
    err_file = open(job_error_path, 'a')
    subprocess.call(shlex.split(work_term), stdout = out_file, stderr = err_file, env=cwd_env, shell=False) #give a file handle; shell=False takes arguments as a list
    err_file.close()
    out_file.close()

    return(0)


def create_result_folder(folder_path=""):
    if (os.path.exists(folder_path)==False):
        os.makedirs(folder_path)
    print("Create working folder: %s" % (folder_path))        
    print("")


def run_ffprofile(program_path="", ext_args_str="", load_path_list=[], save_folder_path="", n_cores=0, cwd_env={}, overwrite_flag=False, job_output_path='', job_error_path=''):

    create_result_folder(save_folder_path)

    args_list=[]

    for load_path in load_path_list:
        save_path = "%s/%s" % (save_folder_path, os.path.basename(load_path))
        work_term = "%s %s %s %s" % (program_path, ext_args_str, load_path, save_path)

        if (overwrite_flag==True or (os.path.exists(load_path) and os.path.getsize(load_path)!=0) ):    
            args_list.append([work_term, job_output_path, job_error_path, cwd_env])
        else:
            print("SKIP: %s" % (work_term))

    if (args_list!=[]):
        spontaneous_process(multiprocess_unit, args_list, n_cores)
    
    print("")


def run_dist_matrix(program_path="", ext_args_str="", load_folder_path="", save_path='', cwd_env={}, overwrite_flag=False, job_output_path='', job_error_path=''):

    create_result_folder(os.dirname(save_path))

    file_list=[]

    for root, dirs, files in os.walk(load_folder_path): #recursive search for (ffprofile) files
        for name in files:
            file_list.append("%s/%s" % (root, name))

    #ext_args_str = "\s".join([])

    work_term = "%s %s %s > %s" % (program_path, ext_args_str, ", ".join(file_list), save_path)

    if (overwrite_flag==True or os.path.exists(save_path)==False or os.path.getsize(save_path)==0):
        multiprocess_unit(work_term, job_output_path, job_error_path, cwd_env)
    else:
        print("SKIP: %s" % (work_term))

    print("")


# https://dendropy.readthedocs.io/en/main/library/phylogeneticdistance.html
# https://dendropy.readthedocs.io/en/main/primer/phylogenetic_distances.html
def dendropy_DistMatrix_to_Tree(load_path="", save_path="", tree_method=""):
    #output: Bool

    pdm = dendropy.PhylogeneticDistanceMatrix()
    tree = ""

    with open(load_path) as read_f: #cannot accept low triangular matrix
        print(read_f.readline()) #skip the first row

        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=read_f
            , is_first_row_column_names=False
            , is_first_column_row_names=True #taxon names enlisted on the first column
            , is_allow_new_taxa=True
            , delimiter="\t" #reason distance matrix should be tab or comma delimited (unsupport delimit_whitespace)
            #, delimit_whitespaces=True
            )
    
    #print(pdm) #dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix object
    #print(pdm.as_data_table) #dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix object
    #print(pdm.taxon_namespace)

    if (tree_method=="nj"):
        tree = pdm.nj_tree()
        #print(tree.as_string("newick"))

    elif (tree_method=="upgma"):
        tree = pdm.upgma_tree()
        #print(tree.as_string("newick"))

    else:
        print("Unsupport: %s" % (tree_method))
        return(False)

    write_f = open(save_path, "w")
    write_f.write(tree+"\n")
    write_f.close()


    print("")


# https://medium.com/geekculture/phylogenetic-trees-implement-in-python-3f9df96c0c32
# https://biopython.org/wiki/Phylo
def bio_DistMatrix_to_Tree(load_path="", save_path="", tree_method=""):
    #output: Bool

    taxon_list=[]
    mat_list=[]
    n_cnt=0

    ## remodel input distance matrix to Bio.Phylo.DistanceMatrix()
    with open(load_path, 'r') as read_f:
        n_item = int(read_f.readline().strip())

        for read_line in read_f:
            n_cnt+=1
            str_list = read_line.strip().split()

            taxon_list.append(str_list[0])
            mat_list.append([float(n_item) for n_item in str_list[1:n_cnt+1]])

    dm = DistanceMatrix(names=taxon_list, matrix=mat_list) #manually create a Bio DistanceMatrix object
    #print(dm)

    tree = DistanceTreeConstructor()

    if (tree_method=="nj"):
        tree = DistanceTreeConstructor().nj(dm)

    elif (tree_method=="upgma"):
        tree = DistanceTreeConstructor().upgma(dm)

    else:
        print("Unsupport: %s" % (tree_method))
        return(False)
    
    Phylo.write(tree, save_path, "newick") #save tree in Newick format

    ''' #using StringIO(). this doesn't work
    str_io = io.StringIO()
    str_io.seek(0)
    save_path = "./save_tree"
    Phylo.write(tree, str_io, "newick") #does not work
    print(str_io.getvalue())
    '''

## short coming of Robinson-Foulds distance
# https://ms609.github.io/TreeDist/articles/Generalized-RF.html
def RF_topology_distance(tree_path_a='', tree_path_b='', tree_format='newick', normalize_flag=False):

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


def argument_str_struct(args_dict={}): #to a single string of arguments

    argument_term=""

    for n_opt in args_dict:

        argument_term+="%s " % (args_dict[n_opt])

    return argument_term.strip()
                      


def read_sym_matrix(load_path=''): #2019-11; can read symmetric and lower-triangular matrix
    n_items = int(pd.read_fwf(load_path, header=None, infer_nrows=1)[0].tolist()[0]) #get number of items = column size

    df = pd.read_csv(
        load_path
        , delim_whitespace=True
        , lineterminator='\n'
        , header=None
        , names=[n_index for n_index in range(0, n_items+1)] #assign temporary column names to the size of n_items
        , index_col=0
        , skiprows=1
        , skip_blank_lines=True
        #, engine="pyarrow"
        ) #, keep_default_na=True) #load file to pd.DataFrame()

    df.index = df.index.astype(str)
    df.columns = df.index #replace column names to df.index

    if (df.isnull().sum(axis=1).values[0]!=0): #transpose and sum
        df = df.fillna(0)
        df+=df.T

    #print(df.columns)
    #print(df.index)
    #print(df.shape)

    return (df) #index list, dataframe


def root_determine(pd_df = pd.DataFrame()):
    ## potential outgroup for rooting purpose. Select the OTU(s) far from the rest.
    rsum = pd_df.sum(axis=1) #so far result the same using mean - std
    root_list = list(pd_df.index[rsum==max(rsum)])
    rsum_mean = statistics.mean(rsum)

    print("\toutgroup/root determine")
    print("Maximum row sum = %.2f, mean row sum: %.2f. suggested root(s) %s" % (max(rsum), rsum_mean, ",".join(root_list)))
    
    return(root_list)


def rooting_tree(newick_str='', reroot_clade_list=[], rooting_method=""): #using dendropy

    # premodification of newick string for convenience
    newick_str = newick_str.replace(" ", "_")

    if (reroot_clade_list!=[]):
        reroot_clade_list = [n_newick.replace("_", " ") for n_newick in reroot_clade_list] #dendropy natively convert "_" within labels to " "

        #print reroot_clade_list
        tree = dendropy.Tree.get_from_string(newick_str, "newick")
        tree.encode_bipartitions()
        #tree.is_rooted = True

        #https://dendropy.readthedocs.io/en/v3.12.1/tutorial/treemanips.html
        #reroot options: edge, node, and outgroup (single)
        if (len(reroot_clade_list)>1):

            mrca = tree.mrca(taxon_labels=reroot_clade_list) #mrca = most recent common ancestor of multiple outgroups
            
            if (rooting_method=="node"):
                tree.reroot_at_node(mrca, update_bipartitions=True)
            
            elif (rooting_method=="edge"):
                half_edge = mrca.edge_length / 2.0
                tree.reroot_at_edge(mrca.edge, length1=half_edge, length2=half_edge, update_bipartitions=True, suppress_unifurcations=False)
                #tree.reroot_at_edge(mrca.edge, update_bipartitions=True)

            else:
                print("unsupported rooting method using outgroups of ", reroot_clade_list)

            ## midpoint rooting unsupported for ultiple outgroups

        else: #pull one outgroup, then halved its edge
            #print(tree.taxon_namespace)
            outgroup_node = tree.find_node_with_taxon_label(reroot_clade_list[0])

            if (rooting_method=="node"):
                tree.to_outgroup_position(outgroup_node, update_bipartitions=True)
            
            elif (rooting_method=="edge"):
                half_edge = outgroup_node.edge_length / 2.0
                tree.reroot_at_edge(outgroup_node.edge, length1=half_edge, length2=half_edge, update_bipartitions=True, suppress_unifurcations=False)
            
            elif (rooting_method=="midpoint"):
                tree.reroot_at_midpoint(update_bipartitions=True) #find the longest edge and then split (not halved)

            else:
                print("unsupported rooting method using an outgroup of ", reroot_clade_list)

        tree.is_rooted = True

        reroot_clade_list = [n_newick.replace(" ", "_") for n_newick in reroot_clade_list] #recover " " to "_"

        newick_str = tree.as_string(schema='newick') #rerooted newick start with [&R]. remove or keep?

    # rollback previous newick string modification
    newick_str = newick_str.replace(" ", "_")
    newick_str = newick_str.replace("'", "") #remove quotations generated during the process

    return(newick_str)



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
            , 'hvs:e:S:m:rR:M:H:Q:'
            , ["ffp=", "ffp_xarg=", "mat=", "mat_xarg=", "rf_dev=", "metric_ubound=", "tree_ext=" ,"tree_xarg=", "tree_method=", "rf_norm", "outgroup=", "root_method="]
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

    for opt, arg in opts:

        if (opt=="-h"):
            show_help()

        elif (opt=="-v"):
            show_version()

        elif (opt=="-s"):
            l_begin=int(arg)

        elif (opt=="-e"):
            l_end=int(arg)

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
            args_dict["-w"]=True

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
    dist_mat_program_path=""
    tree_program_path=""

    rf_distance_list=[]

    cwd_env = os.environ #currento working environment for shell execution

    if (len(args)>0):

        for n_path in args: #can be either files or folders
            load_path = os.path.abspath(n_path)

            #prevent file path interruption in linux environment
            load_path = load_path.replace('(', '\(')
            load_path = load_path.replace(')', '\)')
            load_path = load_path.replace(' ', '\ ')
    
            load_path_list.append(load_path)


        if ((ffp_program_path!="" and os.path.exists(ffp_program_path)) and (dist_program_path!="" and os.path.exists(dist_program_path))):
                
            for feature_length in range(l_begin, l_end+1): #0 base

                if (l_begin==0):
                    save_folder_path = save_work_path

                else:
                    save_folder_path = "%s.%d" % (save_pre_folder_path, feature_length)
                    args_dict["-s"] = "-s %d" % (cy1)


                # initiate paths (save and load)
                fp_save_path = "%s/fp.l.%d" % (save_work_path, feature_length)
                distmat_save_path = "%s/distance_matrix/mat.l.%d" % (save_work_path, feature_length)
                tree_save_path="%s/tree/tree.l.%d" % (save_work_path, feature_length)

                # 1. FFprofile
                run_ffprofile(
                    program_path=ffp_program_path
                    , ext_args_str=args_dict["--ffp_xarg"]
                    , load_path_list=load_path_list
                    , save_folder_path=fp_save_path #folder path
                    , n_cores=n_cores
                    , cwd_env=cwd_env
                    , overwrite_flag=args_dict["-w"]
                    , job_output_path=job_output_path
                    , job_error_path=job_error_path
                    )

                if (l_begin!=0):
                    # 2. dist matrix
                    run_dist_matrix(
                        program_path=dist_mat_program_path
                        , ext_args_str=args_dict["--mat_xarg"]
                        , load_folder_path=fp_save_path
                        , save_path=distmat_save_path
                        , cwd_env=cwd_env
                        , overwrite_flag=args_dict["-w"]
                        , job_output_path=job_output_path
                        , job_error_path=job_error_path
                        )


                    # 3. generate a phylotree
                    if (tree_program_path!="" and os.path.exists(tree_program_path)): #use distance matrix -> tree building program directed by a user.
                        work_term = "%s %s %s %s" % (
                            tree_program_path
                            , args_dict["--tree_xarg"]
                            , distmat_save_path
                            , tree_save_path
                            ).strip()

                        print("Run: %s" % (work_term))

                        multiprocess_unit(work_term, job_output_path, job_error_path, cwd_env, job_output_path, job_error_path)
                    
                    else: #built-in tree recontruction function
                        bio_DistMatrix_to_Tree(
                            load_path=distmat_save_path
                            , save_path=tree_save_path
                            , tree_method=args_dict["--tree_method"]
                            )
                        #dendropy_DistMatrix_to_Tree(load_path="", save_path="", tree_method="")


                    # 4. Robinson-Foulds distance: RF variation to determine to continue search or not
                    #if (feature_length>l_begin): #Rf between consequtive ls (e.g., l-1 and l)
                    previous_tree_save_path="%s/tree/tree.l.%d" % (save_work_path, feature_length-1)
                    # current_tree_save_path = tree_save_path

                    if (os.path.exists(tree_save_path) and os.path.exists(previous_tree_save_path)):
                        #determine the RF distance between trees of l-1 (previous) and l (current), ideally rf_distance==0 is best (reached topological convergence)
                        rf_distance = RF_topology_distance(
                            tree_path_a=previous_tree_save_path
                            , tree_path_b=tree_save_path
                            , tree_format='newick'
                            , normalize_flag=args_dict["--rf_norm"]
                            )
                        
                        rf_distance_list.append(rf_distance)

                        if (abs(rf_distance_list[-2] - rf_distance_list[-1]) > rf_diff_lthreshold): #if RF change is drastic then continue until the different is below the threshold
                            # increment and continue searching for optimum-l
                            l_end+=1

                            ## also check if matrix has the maximum JSD value to prevent searching endless and meaningless
                    else:
                        continue


                    # 5. root tree, get rooted tree inaddition to unroot tree
                    if (reroot_clade_list==[]): #auto determine for rooting
                        reroot_clade_list = root_determine(
                            pd_df = read_sym_matrix(load_path = distmat_save_path))

                    rooted_tree_save_path = "%s.rooted" % (tree_save_path)
                    rooted_tree = ""

                    with open(tree_save_path, 'r') as read_f:
                        rooted_tree = rooting_tree(
                            newick_str=read_f.read()
                            , reroot_clade_list=reroot_clade_list
                            , rooting_method=root_method
                            )

                    with open(rooted_tree_save_path, 'w') as write_f: #save rooted tree
                        write_f.write(rooted_tree+"\n")


        else:
            print("check FF-Profiling program path: %s" % (ffp_program_path))
            print("or check dist matrix program path: %s" % (dist_program_path))

    else:
        print("without input files")

