# python3.12 <=

### README
# You may modifiy and redistribute, at your own risk, to utilize this script not limited to run tutorial.
###

import os, sys, shutil, argparse, io
import subprocess, shlex, time

from time import gmtime, strftime

from multiprocessing import cpu_count, Pool

import dendropy
from dendropy.calculate import treecompare

import pandas as pd #pandas use matlablib as a backend for plotting
import numpy as np

import statistics

#DEFAULT SET
DEFAULT_DEDROPY_TREE_METHOD_LIST=["nj", "upgma"]
MAXIMUM_DISTANCE = 0.999999 #consider floating point precision; practically 1.0

class aux_core:

    def create_folder(
        folder_path=""
        , overwrite_flag=False
        ):

        if (os.path.exists(folder_path)==True and overwrite_flag==True):
            shutil.rmtree(folder_path)

        os.makedirs(folder_path, exist_ok=True)

        print(f"\\/Create working folder: {folder_path} (overwrite={overwrite_flag})")


    def run_multiprocess(
        proc_name
        , args_list=[]
        , n_cores=0
        ):

        result_list=[]

        if (args_list==[]):
            print("\\No work term submitted; pass")

        else:
            
            with Pool(processes=n_cores) as proc_pool: #define pool
                ret_list = proc_pool.starmap_async(
                    func=proc_name
                    , iterable=args_list
                    , chunksize=max(1, int(len(args_list) / n_cores))
                    )

                proc_pool.close()
                proc_pool.join()

                result_list = ret_list.get()


            for index, n_item in enumerate(result_list):
                if (n_item==0): #fail to retrieve or process result
                    print(f"ERROR returning: {args_list[index][0]}")

            if (0 in result_list):
                raise ValueError("#HALT: unable to proceed next steps")

        #return ret_list


    def subprocess_run(
        work_term=''
        , job_output_path=''
        , job_error_path=''
        , std_check_flag=False
        , cwd_env={}
        ):

        try:
            out_file = open(job_output_path, 'w')
            err_file = open(job_error_path, 'a')

            subprocess.run( #advised to use .run() beyond Python 3.3
                shlex.split(work_term)
                , stdout = out_file
                , stderr = err_file
                #, capture_output=True #not compatible with stdout, stderr
                , check=std_check_flag #True wait for stdout or stderr
                , env=cwd_env
                , shell=False
                ) #give a file handle; shell=False takes arguments as a list
            
            err_file.close()
            out_file.close()

            return(1)
        
        except:
            return(0)


    def read_dmatrix(
        load_path=""
        ):

        with open(load_path, 'r') as read_f:
            taxon_cnt = int(read_f.readline()) #skip the first line

            pd_df = pd.read_csv(
                read_f,
                sep=r"\s+", #all white spaces, either spaces or tab
                lineterminator='\n',
                header=None,
                names=[n_index for n_index in range(0, taxon_cnt+1)], #assign temporary column names to the size of n_items
                index_col=0,
                #skiprows=0,
                skip_blank_lines=True,
                low_memory=False
            ) #load file to pd.DataFrame(); missing upper triangular portion of a lower triangular matrix is NaN

        pd_df = pd_df.astype(float)
        pd_df.columns = pd_df.index

        # lower triangular to square (symmetric) matrix; a square matrix is also accepted, using its lower triangular portion
        pd_df = pd_df.where(np.tril(np.ones(pd_df.shape, dtype=bool), k=-1), 0.0) #keep lower triangle, excluding diagonal (is self; 0)
        pd_df = pd_df + pd_df.T.values #mirror lower to upper

        return(pd_df)



class ffp_core: #only generate work_terms to be executed through subprocess
    def profiling(
        program_path=""
        , ext_args_str=""
        , feature_length=1 #explicit feature length (-s) per call
        , load_path_list=[]
        , save_folder_path=""
        , cwd_env={}
        , overwrite_flag=False
        , job_output_path=''
        , job_error_path=''
        , std_check_flag=False
        ):

        args_list=[]

        for load_path in load_path_list:
            save_path = f"{save_folder_path}/{os.path.basename(load_path)}"
            work_term = f"{program_path} -s {feature_length} {ext_args_str} {load_path} {save_path}"

            if (overwrite_flag==False and (os.path.exists(save_path) and os.path.getsize(save_path)!=0)):    
                print(f"SKIP> {save_path} (exists)")
            
            else:
                args_list.append([work_term, job_output_path, job_error_path, std_check_flag, cwd_env])

        return(args_list)


    def pairwise_dmatrix(
        program_path=""
        , ext_args_str=""
        , load_folder_path=""
        , save_path='' #, save_folder_path=""
        , n_cores=0
        , cwd_env={}
        , overwrite_flag=False
        , job_output_path=''
        , job_error_path=''
        , std_check_flag=False
        ):

        args_list=[]
        file_list=[]

        for root, dirs, files in os.walk(load_folder_path): #recursive search for (ffprofile) files
            for name in files:
                file_list.append("%s/%s" % (root, name))

        if (ext_args_str.find("-t")==-1): #if not found (multithreaded)
            ext_args_str+=f" -t {n_cores}" #multithreaded

        if (ext_args_str.find("-T")==-1): #tab delimited; default PHYLIP format truncates taxon names to 9 characters
            ext_args_str+=" -T"

        #work_term = "%s %s %s > %s" % (program_path, ext_args_str, " ".join(file_list), save_path) #used to work with .call()
        work_term = "%s %s %s" % (program_path, ext_args_str, " ".join(file_list).strip())

        if (overwrite_flag==False and (os.path.exists(save_path) and os.path.getsize(save_path)!=0)):    
            print(f"SKIP> {save_path} (exists)")
        
        else:
            args_list.append([work_term, job_output_path, job_error_path, std_check_flag, cwd_env])

        return(args_list)
    

class tree_handle:
    
    # https://dendropy.readthedocs.io/en/main/library/phylogeneticdistance.html
    # https://dendropy.readthedocs.io/en/main/primer/phylogenetic_distances.html
    def dendropy_DistMatrix_to_Tree(
        load_path=""
        , save_path=""
        , tree_method="nj"
        , tree_schema="newick"
        , overwrite_flag=False
        ):

        if (os.path.exists(save_path)==True
            and os.path.getsize(save_path)!=0
            and overwrite_flag==False
            ):

            print(f"SKIP: found an existing non-zero tree file, {save_path}")

        else:   

            taxon_cnt=0
            tree = ""

            pdm = dendropy.PhylogeneticDistanceMatrix()
            str_io = io.StringIO()

            # dendropy cannot accept low triangular matrix; read_dmatrix converts it to a square matrix
            aux_core.read_dmatrix(load_path).to_csv(
                str_io
                , sep="\t"
                , header=False
                , index=True #taxon names on the first column
                )
            str_io.seek(0)

            # require to remove any trailing blank lines at the end
            pdm = dendropy.PhylogeneticDistanceMatrix.from_csv( #considering only the upper right triangle portion
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
                print(f"# Unsupported tree method: {tree_method}")
                return(False)

            if (save_path!=""):
                with open(save_path, "w") as write_f:
                    write_f.write(tree.as_string(tree_schema))


    def external_build_Tree(
        load_path=""
        , save_path=""
        , argument_str=""
        , overwrite_flag=False
        , cwd_env={}
        ):

        work_term = f"{argument_str} {load_path} {save_path}".strip()

        if ((os.path.exists(save_path)==True 
             and os.path.getsize(save_path)!=0)
             and overwrite_flag==False
             ):
            
            print(f"SKIP> {save_path} (exists)")

        else:    
            aux_core.subprocess_run(
                work_term = work_term
                , job_output_path = job_output_path
                , job_error_path = job_error_path
                , std_check_flag=False
                , cwd_env = cwd_env
            )


    def reroot_tree( #simplify the code
        load_path=""
        , save_path=""
        , mrca_clade_list=[]
        , rooting_method=""
        , tree_schema="newick"
        ): #using dendropy

        tree = dendropy.Tree.get_from_path(
            src = load_path
            , schema = tree_schema
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
            outgroup_node = tree.find_node_with_taxon_label(mrca_clade_list[0])

            if (rooting_method=="node"):
                tree.to_outgroup_position(outgroup_node, update_bipartitions=True)
            
            elif (rooting_method=="edge"):
                half_edge = outgroup_node.edge_length / 2.0
                tree.reroot_at_edge(outgroup_node.edge, length1=half_edge, length2=half_edge, update_bipartitions=True, suppress_unifurcations=False)

            else:
                raise ValueError(f"unsupported rooting method: {rooting_method} using {mrca_clade_list}")

        tree.is_rooted = True

        with open(save_path, "w") as write_f:
            write_f.write(tree.as_string(tree_schema))


    ## short coming of Robinson-Foulds distance
    # https://ms609.github.io/TreeDist/articles/Generalized-RF.html
    def RF_distance(
        tree_path_a=''
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
            #normalized robinson-fould distance; 2n - 6 is the maximum difference if unrooted and binary trees. n is a number of tips or taxons
            return float(rf_distance) / float(len(tree_a.taxon_namespace) + len(tree_b.taxon_namespace) - 6)

        else:
            return rf_distance


class do_plot: #combine RF distance and upperbound in one plot
    def search_optimum_range(
        pd_df = pd.DataFrame()
        , x_axis=""
        , y_axis=""
        , is_y_axis_normalized=False
        , save_path=""
        , at_upperbound=0
        ):

        ax = pd_df.plot.line(
            x=x_axis
            , y=y_axis
            , xticks = list(pd_df[x_axis]) #fix to integer
            , title = "Optimum l search space"
            , xlabel = f"{x_axis} (l versus l+1)"
            , ylabel = f"{y_axis} (normalized)" if is_y_axis_normalized==True else y_axis
            , grid = True
            , legend = False
        )

        if (at_upperbound > 0): #0 means not found
            ax.axvline( ## show upperbound in a vertical line
                x=at_upperbound
                , color='r'
                , linestyle='--'
                )

        ax.get_figure().savefig(
            save_path
            , bbox_inches="tight"
            ) #save to .png



def root_determine(pd_df = pd.DataFrame()): #numerically determine outgroup from a distance matrix; the taxon most distant from the rest
    rsum = pd_df.sum(axis=1) #so far result the same using mean - std
    root_list = list(pd_df.index[rsum==max(rsum)])
    rsum_mean = statistics.mean(rsum)

    print("\toutgroup/root determine")
    print("Maximum row sum = %.2f, mean row sum: %.2f. suggested root(s) %s" % (max(rsum), rsum_mean, ",".join(root_list)))
    
    return(root_list)


def build_arg_parser(): # replaces previous getopt.getopt() and show_help()
    parser = argparse.ArgumentParser(
        description='Run feature frequency profiling (FFP) up to calculating distance matrix. Run "profile", then "distance" with the FF-Profile folders, [-S]/fp.[l], then "optimum" with the distance matrices, [-S]/distance_matrix/matrix.[distance type].[l]'
        )
    parser.add_argument('-v', '--version', action='store_true', help='show version and citation information')

    # feature length range, for profile
    range_parser = argparse.ArgumentParser(add_help=False)
    feature_range = range_parser.add_argument_group("Feature length iteration from [-s] to [-e]")
    feature_range.add_argument('-s', dest='l_begin', type=int, default=1, metavar='INT', help='from feature length (default=1)')
    feature_range.add_argument('-e', dest='l_end', type=int, default=1, metavar='INT', help='to feature length (default=1)')

    # options shared by all subcommands
    common_parser = argparse.ArgumentParser(add_help=False)
    common = common_parser.add_argument_group("Process behavior")
    common.add_argument('-S', dest='work_dir_path', default=os.getcwd(), metavar='PATH', help='a whole process and all output are constructed in the designated folder')
    common.add_argument('-w', '--overwrite', dest='overwrite_flag', action='store_true', help='overwrite existing output files of the same path; without it, existing (non-empty) output files are skipped')
    common.add_argument('-t', '--thread', dest='n_cores', type=int, default=max(1, cpu_count() -1), metavar='INT', help='a number of threads or cores to use')

    subparsers = parser.add_subparsers(dest='command', metavar='{profile,distance,optimum}')

    ## 1. Construct FF-Profile
    profile_parser = subparsers.add_parser('profile', parents=[range_parser, common_parser]
        , help='1. Construct FF-Profiles; saved to [-S]/fp.[l]'
        , description='1. Construct FF-Profiles of input files for each feature length; saved to [-S]/fp.[l]')
    profile = profile_parser.add_argument_group("1. Construct FF-Profile")
    profile.add_argument('--ppath', required=True, metavar='PATH', help='the path to FF-Profiling program')
    profile.add_argument('--parg', default='', metavar='STR', help='arguments for FF-Profiling program except input and output')
    profile_parser.add_argument('input_paths', nargs='+', metavar='path_to_input_files', help='input sequence files')

    ## 2. Construct distance matrix
    dist_parser = subparsers.add_parser('distance', parents=[common_parser]
        , help='2. Construct (pairwise) distance matrix from FF-Profile folders (fp.[l]); saved to [-S]/distance_matrix/matrix.[distance type].[l]'
        , description='2. Construct (pairwise) distance matrix from FF-Profile folders made by "profile", sorted by feature length [l] (the numeric folder name suffix) in ascending order; '
            'saved to [-S]/distance_matrix/matrix.[distance type].[l]. Distance type is taken from -d in --parg (default=jsdiv)')
    dist_parser.add_argument('input_paths', nargs='+', metavar='path_to_profile_folders'
        , help='FF-Profile folders, e.g., [-S]/fp.*')
    dist_parser.add_argument('-f', '--force', dest='search_enforce_flag', action='store_true'
        , help='force to search all given feature lengths beyond upperbound. Upperbound is a point where the first maximum distance (between ingroup taxons) found')
    dist_parser.add_argument('--outgroup', default='', metavar='STR'
        , help='designate outgroup taxon; multiple outgroups separated by comma. Specified outgroup(s) are excluded when determining upperbound (of ingroups)')

    dmat = dist_parser.add_argument_group("2. Construct (pairwise) Distance Matrix")
    dmat.add_argument('--ppath', required=True, metavar='PATH', help='the path to dist_matrix program')
    dmat.add_argument('--parg', default='', metavar='STR', help='arguments for dist_matrix program except input and output; -t [thread] and -T (tab delimited) are added if missing')

    ## 3-4. Construct tree, RF distance, summary and plot
    opt_parser = subparsers.add_parser('optimum', parents=[common_parser]
        , help='3-4. Construct tree from distance matrices (matrix.[distance type].[l]), then RF distance, summary and plot'
        , description='3-4. Construct (rooted) tree from distance matrices made by "distance", sorted by feature length [l] (the numeric file name suffix) in ascending order; '
            'calculate RF distance between trees of adjacent ls, then summary and plot. Trees are saved to [-S]/tree/tree.[distance type].[l]')
    opt_parser.add_argument('--outgroup', default='', metavar='STR'
        , help='designate outgroup taxon; multiple outgroups separated by comma. Specified outgroup(s) will be used for rooting trees and determining upperbound (of ingroups)')
    opt_parser.add_argument('input_paths', nargs='+', metavar='path_to_distance_matrices'
        , help='distance matrix files, e.g., [-S]/distance_matrix/matrix.jsdiv.*; lower triangular or square matrix')

    tree = opt_parser.add_argument_group("3. Construct Tree")
    tree.add_argument('--ppath', default='', metavar='PATH', help='the path to (external) tree build program of choice. Without --ppath will use dendropy for tree construction using --method')
    tree.add_argument('--parg', default='', metavar='STR', help='arguments for (external) tree build program except input and output. Will call external tree build program: [ppath] [parg] [input] [output], in this order')
    tree.add_argument('--method', default='nj', type=str.lower, choices=DEFAULT_DEDROPY_TREE_METHOD_LIST, help='dendropy tree construction method (default=nj)')
    tree.add_argument('--root_method', default='node', type=str.lower, choices=['node', 'edge'], help='rooting by "node" or "edge" (default=node)')
    tree.add_argument('--schema', default='newick', type=str.lower, choices=['newick', 'nexus', 'nexml'], help='tree format types (default=newick)')

    summary = opt_parser.add_argument_group("4. Summary and Plot")
    summary.add_argument('--rf_norm', action='store_true', help='normalize RF distance; recommended to use (True); default is False')

    return(parser)


def join_parg_values(argv_list):
    # argparse rejects a value starting with '-' (e.g. --parg "-n"), which getopt accepted.
    # Join pass-through arguments as --parg=value to keep the same command line working.
    parg_opt_set = {"--parg"}
    joined_list = []
    n_index = 0

    while (n_index < len(argv_list)):
        if (argv_list[n_index] in parg_opt_set and n_index+1 < len(argv_list)):
            joined_list.append(f"{argv_list[n_index]}={argv_list[n_index+1]}")
            n_index += 2

        else:
            joined_list.append(argv_list[n_index])
            n_index += 1

    return(joined_list)



def sort_by_feature_length(path_list=[], expect_format=""):
    # input files or folders named [prefix].[l]; return [(l, absolute path), ...] sorted by l (numeric suffix) in ascending order
    load_list=[]

    for n_path in path_list:
        n_path = os.path.normpath(n_path) #remove a trailing slash of a folder path
        l_suffix = os.path.basename(n_path).rsplit(".", 1)[-1]

        if (l_suffix.isdigit()==False):
            raise ValueError(f"# Feature length (numeric suffix) not found in the name: {n_path}; expect {expect_format}")

        load_list.append((int(l_suffix), os.path.abspath(n_path)))

    load_list.sort() #by feature length

    return(load_list)


def show_version():
    print('FFP-Wrapper, version 2.0 (since Jan-2025)')
    print("\tDetails and citation information at https://github.com/jaejinchoi/FFP/tree/master")
    print("\tTutorial available at https://github.com/jaejinchoi/FFP/tree/master/example")
    
    sys.exit()


if __name__=='__main__':

    cwd_env = None #os.environ #deliver current working environment flags; not necessary; raise error for pool.starmap

    parser = build_arg_parser()
    opts = parser.parse_args(join_parg_values(sys.argv[1:]))

    if (opts.version==True):
        show_version()

    if (opts.command==None):
        parser.print_help()
        sys.exit()

    # necessary parameters
    l_begin=getattr(opts, "l_begin", 1) #1 as minimum; optimum takes ls from input file names
    l_end=getattr(opts, "l_end", 1) #default=1
    overwrite_flag=opts.overwrite_flag
    work_dir_path=os.path.abspath(opts.work_dir_path)
    n_cores = max(1, opts.n_cores) #least 1

    args_dict={ # arguments with the default values
        "--ffp_ppath": "" #ffprofile program path
        , "--ffp_parg": "" #arguments for ffprofile program

        , "--dmat_ppath": "" #dist matrix program path
        , "--dmat_parg": "" #arguments for dist matrix program

        , "--tree_ppath": "" #(external) tree build program path
        , "--tree_parg": "" #arguments for tree build program
        , "--tree_method": "nj"
        , "--tree_schema": "newick"
        , "--rf_norm": False
        , "--root_method": "node"

    }

    if (opts.command=="profile"):
        args_dict["--ffp_ppath"] = os.path.abspath(opts.ppath) #ffprofile program path
        args_dict["--ffp_parg"] = str(opts.parg).strip() #external arguments

    elif (opts.command=="distance"):
        search_enforce_flag=opts.search_enforce_flag
        outgroup_list = str(opts.outgroup).split(",") if opts.outgroup!="" else [] #exclude outgroup when determining upperbound

        args_dict["--dmat_ppath"] = os.path.abspath(opts.ppath) #dist matrix program path
        args_dict["--dmat_parg"] = str(opts.parg).strip() #external arguments

        # distance type for output file names; value of -d (or --distance) in --dmat_parg, else the program default
        dmat_type = "jsdiv"
        dmat_parg_tokens = shlex.split(args_dict["--dmat_parg"])
        for n_index, n_token in enumerate(dmat_parg_tokens[:-1]):
            if (n_token in ["-d", "--distance"]):
                dmat_type = dmat_parg_tokens[n_index+1]

    elif (opts.command=="optimum"):
        outgroup_list = str(opts.outgroup).split(",") if opts.outgroup!="" else [] #designate outgroup for rooting purpose

        args_dict["--tree_ppath"] = os.path.abspath(opts.ppath) if opts.ppath!="" else "" #(external) tree build program path
        args_dict["--tree_parg"] = str(opts.parg).strip() #external arguments
        args_dict["--tree_method"] = opts.method
        args_dict["--tree_schema"] = opts.schema
        args_dict["--root_method"] = opts.root_method
        args_dict["--rf_norm"] = opts.rf_norm #normalize RF distance


    #create working directory if missing
    os.makedirs(work_dir_path, exist_ok=True) #exist_ok=True to avoid error if the folder already exists

    #pre-check
    if (l_begin>0 and l_end<l_begin): #single l run
        l_end = l_begin

    elif (l_begin<0):
        print("# l_start and l_end should be larger than 0, and l_start <= l_end")
        sys.exit()

    # record subprocess output and error
    current_time_str = strftime("%m.%d.%Y", gmtime())
    job_output_path="%s/%s_out.txt" % (work_dir_path, current_time_str)
    job_error_path="%s/%s_error.txt" % (work_dir_path, current_time_str) #can be avoid by using qsub -j y


    if (opts.command=="profile"):

        load_path_list=[] # of FASTA formatted sequence files

        for n_path in opts.input_paths: #can be either files or folders
            load_path = os.path.abspath(n_path)

            #prevent file path interruption in linux environment
            load_path = load_path.replace('(', '\\(')
            load_path = load_path.replace(')', '\\)')
            load_path = load_path.replace(' ', '\\ ')
            load_path_list.append(load_path)

        # 0. Check available programs
        if (os.path.exists(args_dict["--ffp_ppath"])==False):
            print(f"Core program {args_dict["--ffp_ppath"]} not found")
            raise ValueError(f"# Please check if it is compiled and executable")

        if ("-s" in shlex.split(args_dict["--ffp_parg"])):
            print(f'! Warning: -s in --parg ({args_dict["--ffp_parg"]}) overrides the feature length from -s/-e of the wrapper')


        # check iteration time (for each feature length)

        # The main iteration per feature length
        for feature_length in range(l_begin, l_end+1):
            iter_begin_at = time.perf_counter() #reset per feature length
            ## save paths
            output_path_ffp_folder = f"{work_dir_path}/fp.{feature_length}"

            print(f"---For Feature length: {feature_length}---")


            ## 1. Construct Feature Frequency Profiles, takes the argument --ffp_xarg
            print(f"# 1. Construct Feature Frequency Profiles")

            aux_core.create_folder(
                folder_path = output_path_ffp_folder
                , overwrite_flag=False #keep the folder; -w overwrites individual profile files
            )

            args_list = ffp_core.profiling(
                program_path=args_dict["--ffp_ppath"]
                , ext_args_str=args_dict["--ffp_parg"]
                , feature_length=feature_length #from -s/-e loop
                , load_path_list=load_path_list
                , save_folder_path=output_path_ffp_folder
                , cwd_env=cwd_env
                , overwrite_flag=overwrite_flag #replace old saved files
                , job_output_path=job_output_path
                , job_error_path=job_error_path
                , std_check_flag=False
                )

            aux_core.run_multiprocess(
                aux_core.subprocess_run
                , args_list=args_list
                , n_cores=n_cores
                )

            print(f"---Done\n")


            ## print ellapsed time per cycle in seconds
            print(f"*-- The execute cycle for l={feature_length}, of {len(load_path_list)} taxons, took approximately {(time.perf_counter() - iter_begin_at):.2f} seconds\n")

        print(f"\n---Finish---\n")


    elif (opts.command=="distance"):

        # 0. Check available programs
        if (os.path.exists(args_dict["--dmat_ppath"])==False):
            print(f"Core program {args_dict["--dmat_ppath"]} not found")
            raise ValueError(f"# Please check if it is compiled and executable")


        # input FF-Profile folders, fp.[l]; sort by l (numeric suffix) in ascending order
        ffp_folder_load_list = sort_by_feature_length(opts.input_paths, expect_format="fp.[l]")

        for feature_length, output_path_ffp_folder in ffp_folder_load_list:
            if (os.path.isdir(output_path_ffp_folder)==False):
                raise ValueError(f"# FF-Profile folder not found: {output_path_ffp_folder}")

        # check iteration time (for each feature length); calculating distance matrix

        # The main iteration per feature length
        for feature_length, output_path_ffp_folder in ffp_folder_load_list:
            iter_begin_at = time.perf_counter() #reset per feature length
            ## save paths
            output_path_dmatrix = f"{work_dir_path}/distance_matrix/matrix.{dmat_type}.{feature_length}"

            print(f"---For Feature length: {feature_length}---")


            ## 2. Calculate pairwise distance matrix
            print(f"# 2. Construct pairwise distance matrix")

            aux_core.create_folder(
                folder_path = os.path.dirname(output_path_dmatrix)
                , overwrite_flag=False #shared by all ls; keep the folder, -w overwrites individual files
            )

            args_list = ffp_core.pairwise_dmatrix(
                program_path=args_dict["--dmat_ppath"]
                , ext_args_str=args_dict["--dmat_parg"]
                , load_folder_path=output_path_ffp_folder
                , save_path=output_path_dmatrix
                , n_cores=n_cores #append as -t (thread) option to ext_args_str if not found
                , cwd_env=cwd_env
                , overwrite_flag=overwrite_flag #replace old saved files
                , job_output_path=output_path_dmatrix #job_output_path
                , job_error_path=job_error_path
                , std_check_flag=True
                )

            #print(args_list[0][0])

            #'''
            aux_core.run_multiprocess(
                aux_core.subprocess_run
                , args_list=args_list
                , n_cores=1 #multithreaded in the called program; iterate per [l]
                )
            #'''
            print(f"---Done\n")


            ## 2.1 Read distance matrix
            pd_df_dmatrix = aux_core.read_dmatrix(output_path_dmatrix)

            ### to determine upperbound: when the first maximum distance (between ingroup taxons) appears
            ### consider only an upper triangular (matrix) values, excluding diagonal distance (is self; 0)
            pd_df_dmatrix = pd_df_dmatrix.iloc[~pd_df_dmatrix.index.isin(outgroup_list), ~pd_df_dmatrix.columns.isin(outgroup_list)] #exclude outgroup-paired distances
            pd_df_dmatrix = pd_df_dmatrix.where( #mask in the shape after excluding outgroups
                np.triu(np.ones(pd_df_dmatrix.shape)
                        , k=1 #excluding diagonal distance (is self; 0)
                        ).astype(bool) #replace lower triangular portion to NA
                        ).stack()

            maximum_distance_count = pd_df_dmatrix[pd_df_dmatrix > MAXIMUM_DISTANCE].count()


            ## print ellapsed time per cycle in seconds
            print(f"*-- The execute cycle for l={feature_length} took approximately {(time.perf_counter() - iter_begin_at):.2f} seconds\n")


            ## upperbound check point. Stop unless instructed to enforce to search remaining feature lengths
            if (search_enforce_flag==False and maximum_distance_count>0):
                print(f'! Upperbound found at l={feature_length}; stop search beyond')
                break #skip searching in longer feature lengths

        print(f"\n---Finish---\n")


    elif (opts.command=="optimum"):

        # to build summary dataframe
        feature_length_list=[]
        matrix_path_list=[]
        tree_path_list=[]
        maximum_distance_count_list=[]


        # input distance matrices, matrix.[distance type].[l]; sort by l (numeric suffix) in ascending order
        dmatrix_load_list = sort_by_feature_length(opts.input_paths, expect_format="matrix.[distance type].[l]")

        # check iteration time (for each feature length); from distance matrix to building a tree

        # The main iteration per feature length
        for feature_length, output_path_dmatrix in dmatrix_load_list:
            iter_begin_at = time.perf_counter() #reset per feature length
            ## save paths; tree.[distance type].[l]
            output_path_tree = f"{work_dir_path}/tree/tree.{os.path.basename(output_path_dmatrix).removeprefix('matrix.')}"

            feature_length_list.append(feature_length)
            matrix_path_list.append(output_path_dmatrix)
            tree_path_list.append(output_path_tree)

            print(f"---For Feature length: {feature_length}---")
            ## 2.1 Read distance matrix
            pd_df_dmatrix = aux_core.read_dmatrix(output_path_dmatrix)

            ### to determine upperbound: when the first maximum distance (between ingroup taxons) appears
            ### consider only an upper triangular (matrix) values, excluding diagonal distance (is self; 0)
            pd_df_dmatrix = pd_df_dmatrix.iloc[~pd_df_dmatrix.index.isin(outgroup_list), ~pd_df_dmatrix.columns.isin(outgroup_list)] #exclude outgroup-paired distances
            pd_df_dmatrix = pd_df_dmatrix.where( #mask in the shape after excluding outgroups
                np.triu(np.ones(pd_df_dmatrix.shape)
                        , k=1 #excluding diagonal distance (is self; 0)
                        ).astype(bool) #replace lower triangular portion to NA
                        ).stack()

            maximum_distance_count = pd_df_dmatrix[pd_df_dmatrix > MAXIMUM_DISTANCE].count()
            maximum_distance_count_list.append(maximum_distance_count)

            ## 3. Build a tree from a distance matrix, takes the argument --dmatrix_xarg
            print(f"# 3. Construct tree")

            aux_core.create_folder(
                folder_path = os.path.dirname(output_path_tree)
                , overwrite_flag=False #shared by all ls; keep the folder, -w overwrites individual files
            )

            if (args_dict["--tree_method"] in DEFAULT_DEDROPY_TREE_METHOD_LIST and args_dict["--tree_ppath"]==""):
                ### build a tree using dendropy
                print(f'* Dendropy --tree_method={args_dict["--tree_method"]}, --tree_schema={args_dict["--tree_schema"]}')

                tree_handle.dendropy_DistMatrix_to_Tree(
                    load_path=output_path_dmatrix
                    , save_path=output_path_tree
                    , tree_method=args_dict["--tree_method"]
                    , tree_schema=args_dict["--tree_schema"]
                    , overwrite_flag=overwrite_flag
                    )

            elif (args_dict["--tree_ppath"]!="" and os.path.exists(args_dict["--tree_ppath"])==True):
                ### use external tree building program if specified
                print(f'* Dendropy --tree_method={args_dict["--tree_method"]}, --tree_schema={args_dict["--tree_schema"]}')

                tree_handle.external_build_Tree(
                    load_path=output_path_dmatrix
                    , save_path=output_path_tree
                    , argument_str=f"{args_dict["--tree_ppath"]} {args_dict['--tree_parg'].strip()} "
                    , overwrite_flag=overwrite_flag #replace old saved files
                    , cwd_env=cwd_env
                    )

            else:
                print(f"# Dendropy tree method: {args_dict['--tree_method']} not supported. Select among {DEFAULT_DEDROPY_TREE_METHOD_LIST}")
                print(f"# OR check --tree_ppath={args_dict["--tree_ppath"]}")
                sys.exit()


            ## 3.a Reroot the tree if instructed, takes the argument --outgroup
            if (outgroup_list!=[]):
                print(f'* Rooting tree. Set outgroup(s) = {outgroup_list}')

                tree_handle().reroot_tree(
                    tree_load_path=output_path_tree
                    , tree_save_path=output_path_tree #overwrite the original tree
                    , mrca_clade_list=outgroup_list
                    , rooting_method=args_dict["--root_method"]
                    , tree_schema=args_dict["--tree_schema"]
                    )
            # else, maintain unrooted tree
            print(f"---Done\n")


            ## print ellapsed time per cycle in seconds
            print(f"*-- The execute cycle for l={feature_length} took approximately {(time.perf_counter() - iter_begin_at):.2f} seconds\n")

        # end of feature length iteration
        print(f"\n---Conclude---")

        # 4. Calculate RF distance (tree of l and l+1); create summary; save plot
        print(f"# 4. Calculate RF distance, create summary, and save plot")
        ## 4.a. Robinson-Foulds distance between trees of two adjacent feature lengths

        print(f"-.a Calculate Robinson-Foulds distance between trees of l and l+1")

        rf_distance_list=[]

        for n_index in range(0, len(feature_length_list)-1):
            rf_distance = tree_handle.RF_distance(
                tree_path_a=tree_path_list[n_index]
                , tree_path_b=tree_path_list[n_index+1]
                , tree_format=args_dict["--tree_schema"] #same format as saved trees
                , normalize_flag=args_dict["--rf_norm"] #recommended True; default is False
                )

            rf_distance_list.append(float(rf_distance))

        rf_distance_list.append(np.nan) #last item; np.nan accepted in numeric calculation, such as min()
        #print(tree_path_list)

        ## 4.b. Construct a summary dataframe
        print(f"-.b Construct a summary dataframe")
        summary_pd_df = pd.DataFrame()
        summary_save_path = f"{work_dir_path}/summary.tsv"

        summary_pd_df["Feature length"] = feature_length_list
        summary_pd_df["matrix_save_path"] = matrix_path_list
        summary_pd_df["tree_save_path"] = tree_path_list
        summary_pd_df["Maximum distance count (without outgroups)"] = maximum_distance_count_list
        summary_pd_df["RF distance"] = rf_distance_list

        ## ** Search upperbound; upperbound based on distance matrix is determined below
        at_upperbound = summary_pd_df.loc[summary_pd_df["Maximum distance count (without outgroups)"]>0, "Feature length"].iloc[0] if sum(maximum_distance_count_list) > 0 else 0 #zero means not found in search range
        at_min_rf_distance = summary_pd_df.loc[summary_pd_df["RF distance"]==min(rf_distance_list), "Feature length"].iloc[0]

        with open(summary_save_path, 'w') as write_f:
            write_f.write(f"# Upperbound at l={at_upperbound}\n")
            write_f.write(f"# Minimum RF distance at l={at_min_rf_distance}\n")

            summary_pd_df.to_csv(
                write_f
                , header=True
                , index=False
                , sep="\t"
                )
        print(f"--saved: {summary_save_path}")

        ## 4.c. Create plot for optimum range search
        print(f"-.c Construct a plot for optimum range search")
        optimum_plot_save_path = f"{work_dir_path}/search_optimum.png"

        do_plot.search_optimum_range(
            pd_df=summary_pd_df
            , x_axis = "Feature length"
            , y_axis = "RF distance"
            , is_y_axis_normalized = args_dict["--rf_norm"]
            , save_path=optimum_plot_save_path
            , at_upperbound=at_upperbound
            )

        print(f"--saved: {optimum_plot_save_path}")

        print(f"\n---Finish---\n")
