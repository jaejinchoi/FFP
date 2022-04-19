# ---
# title: "Visually_determine_k-mer_optimum"
# author: "JaeJin Choi"
# date: "9/23/2020"
# ---

require(getopt)
require(dplyr)
require(stats)
require(ape) #nj and bionj
require(phangorn) #upgma, RF

require(parallel) #for multi-core
# require(pcaPP) #https://rdrr.io/cran/pcaPP/man/cor.fk.html, for fast Kendall-tau estimation
require(ccaPP) #c++ implemented calculation; corKendall corPearson, corSpearman; https://search.r-project.org/CRAN/refmans/ccaPP/html/corFunctions.html
# require(cli) #load the latest version >=3.1.0

require(permute)

require(optparse) #R optget in python-style
opt_list = list(
  make_option(c("-a", "--all"), type="character", action="store_true", default=F, help="Check and run all projects enlisted in project_load_list (.ods or .txt)")
  , make_option(c("-p", "--path"), type="character", action="store", default="", help="Set working directory")
  , make_option(c("-m", "--matrixpre"), type="character", action="store", default="mat.", help="Set matrix file prefix")
  , make_option(c("-t", "--treepre"), type="character", action="store", default="tree.", help="Set tree file prefix")
  , make_option(c("-o", "--outgroup"), type="character", action="store", default="none", help="Set outgroup(s), delimit with comma','")
  # , make_option(c("-p", "--path"), type="character", action="store", default="", help="Set specific project working directory")
  
  # make_option(c("-m", "--maintitle"), type="character", action="store", default="", help="Set plot's main title"),
  # make_option(c("-s", "--subtitle"), type="character", action="store", default="", help="Set plot's subtitle"),
  # make_option(c("-c", "--caption"), type="character", action="store", default="", help="Set plot's caption/description underneath figure"),
  # make_option(c("-S", "--savefile_path"), type="character", action="store", default="", help="Set savefile name"),
  # make_option(c("-k", "--kit"), type="integer", action="store", default=0, help="k=1 (RF plot; smoothing), k=2 (distanec distribution), k=3 (line graph)"),
  # make_option(c("-e", "--opttemplate"), type="character", action="store", default="", help="value template")
)


df_shuffle <- function(df) #shuffle distance matrix order to stress(test) topology variation
{
  m_df <- as.matrix(df) #temporary convert a.dist() to matrix(), full square distance matrix
  shuffled_item_names <- shuffle(colnames(m_df))
  
  # shuffled_df <- as.dist(m_df[shuffled_item_names, shuffled_item_names])
  return (as.dist(m_df[shuffled_item_names, shuffled_item_names]))
}


matrix_iorder_shuffle <- function(dist_mat) #random shuffle
{
  shuffled_names <- shuffle(colnames(dist_mat))
  return(dist_mat[shuffled_names, shuffled_names])
  
}

matrix_iorder_sort <- function(dist_mat)
{
  dist_mat <- dist_mat[order(rownames(dist_mat)), order(colnames(dist_mat))] #sort row
  return(dist_mat)
  
}


subset_dist_mat <- function(dist_mat, exclude_list)
{
  edist_mat<- dist_mat[
    !(row.names(dist_mat) %in% exclude_list) #row subset
    , !(colnames(dist_mat) %in% exclude_list) #col subset
    ]
  
  return(edist_mat)
}


#' @title dtree_mat_corr
#' @Description Calculate correlation-r of input distance matrix and the tree-derived distance matrix
#' @param input_dist_mat distance_matrix
#' @param phylo_tree phylo tree object
#' @param corr_type specify correlation-r method (pearson, spearman, or kendall)
#' @return correlation-r
#' @examples read_to_dist_df(input_dist_mat, phylo_tree_obj, corr_type="")


#' @title read_to_dist_df
#' @Description Read distance matrix file and output dist object
#' @param load_path path_to_distance_matrices (string)
#' @return as.dist object
#' @examples read_to_dist_df(load_path)
read_matrix_to_df <- function(load_path)
{
  tmp_colnames <-  c(seq(1:(length(readLines(load_path))))) #temporarly assign to fit a number of columns when filled in read.table
  
  input_table <- read.table(file = load_path
                            , sep = "" # all white spaces (reduce)
                            , header= F
                            , skip = 1 #skip the first row (a number of OTUS)
                            , as.is = T #recycling rule to guess variable class
                            , col.names = tmp_colnames
                            , fill = T #able to read triangular matrix
                            
  )
  
  input_df <- as.data.frame(input_table[,-1]) #exclude the first column which is a list of items
  colnames(input_df) <- input_table$X1 #place item names (OTU) as column names
  row.names(input_df) <- input_table$X1 #place item names OTU) as row names
  
  if (rowSums(is.na(input_df[1,]))!=0) #when triangular matrix given, it is necessary to make it symetric before sort by item name order
  {
    input_df[is.na(input_df)] <- 0 #set all NAs to 0
    input_df <- input_df + t(input_df) #add lower triangular + upper triangular
  }
  
  dist_df <- as.dist(dist_matrix_sort(input_df))

  return (dist_df)
  
}


read_dist_mat <- function(load_path)
{
  tmp_colnames <-  c(seq(1:(length(readLines(load_path))))) #temporarly assign to fit a number of columns when filled in read.table
  
  dist_table <- read.table(
    file = load_path
    , sep = "" # all white spaces (reduce)
    , header= F
    , skip = 1 #skip the first row (a number of OTUS)
    , as.is = T #recycling rule to guess variable class
    , col.names = tmp_colnames # to align
    , fill = T #able to read triangular matrix
  )
  
  ## reform to square, then assign names
  dist_df <- as.data.frame(dist_table[,-1]) #exclude the first line, which is a number of OTUs
  # print(dist_df)
  colnames(dist_df) <- dist_table$X1 #place OTU names to column names
  row.names(dist_df) <- dist_table$X1 #place OTU names to row names
  
  if (rowSums(is.na(dist_df[1,]))!=0) #when triangular matrix given, it is necessary to make it symetric before sort by item name order
  {
    dist_df[is.na(dist_df)] <- 0 #set all NAs to 0s before merge
    dist_df <- dist_df + t(dist_df) # lower triangle + upper triangle
  }
  
  # return(as.dist(dist_df)) # return ape::dist() object, say mdist
  return(as.matrix(dist_df)) # return as a matrix
  
}


build_phylo_tree <- function(l_dmat, tree_save_path, tree_algorithm, shuffle_flag, overwrite_flag)
{
  output_tree <- ""
  
  if (shuffle_flag==TRUE)
  {
    l_rdmat <- as.dist(matrix_iorder_shuffle(as.matrix(l_dmat)))
      
  } else
  {
    l_rdmat <- as.dist(l_dmat)
  }
  
  if (tree_algorithm=='bionj')
  {
    output_tree <- ape::bionj(l_rdmat)
    
  } else if (tree_algorithm=='nj')
  {
    output_tree <- ape::nj(l_rdmat)
    
  } else if (tree_algorithm=='upgma')
  {
    output_tree <- phangorn::upgma(l_rdmat)
    
  }
  
  if (!file.exists(tree_save_path) || overwrite_flag==T)
  {
    write.tree(output_tree, file = tree_save_path)  
  }
  
  return (output_tree)
}


## determine (adjacent) topological distance
adjacent_topology_distance <- function(
  pri_l
  , lat_l
  , pri_tree
  , lat_tree
  , method_type
  , adjacent_flag
)
{
  output_value <- NA
  
  # pri_tree <- read.tree(file = pri_tree_path, keep.multi = F)
  # lat_tree <- read.tree(file = lat_tree_path, keep.multi = F)
  
  if (adjacent_flag==T && lat_l - pri_l !=1)
  {
    return(output_value)
    
  } else if (method_type=="topo_rf") #Robinson-Foulds distance
  {
    output_value <- phangorn::RF.dist(pri_tree, lat_tree, check.labels=TRUE, normalize=F)
    
  } else if (method_type=="topo_spr") #Subtree pruning and regrafting distance (usually slower than calculating RF)
  {
    output_value <- phangorn::SPR.dist(pri_tree, lat_tree, check.labels=TRUE)
  } else
  {
    stop(paste0("Unsupported topology distance requested"))
  }
}


## determine (adjacent) matrix correlation coefficient
adjacent_mat_corr <- function(
  pri_l
  , lat_l
  , pri_dmat
  , lat_dmat
  , method_type
  , adjacent_flag
  )
{
    output_value <- NA
    
    if (adjacent_flag==T && lat_l - pri_l !=1)
    {
      return(output_value)
      
    } else if (comp_method_type=="pearson")
    {
      output_value <- ccaPP::corPearson(as.vector(pri_dmat), as.vector(lat_dmat))
      
    } else if (comp_method_type=="spearman")
    {
      output_value <- ccaPP::corSpearman(as.vector(pri_dmat), as.vector(lat_dmat), consistent = FALSE)
      
    } else if (comp_method_type=="kendall")
    {
      output_value <- ccaPP::corKendall(as.vector(pri_dmat), as.vector(lat_dmat), consistent = FALSE)
      
    } else
    {
      stop(paste0("Unsupported correlation method (adjacent_mat_corr) requested"))
    }
    
}
  
## determine correlation coefficient between an original and the tree-derived pariwise distance matrices
mat_dtree_corr <- function(
  l_dmat
  , l_tree
  , method_type
)
{
  output_value <- NA
  
  ## https://rdrr.io/cran/ape/man/cophenetic.phylo.html; calculate tree nodes distances
  ## calculate pariwise distance matrix from given tree
  dtree_dmat <- as.dist(matrix_isort_sort(cophenetic.phylo(l_tree)))
  ## http://www2.uaem.mx/r-mirror/web/packages/ccaPP/ccaPP.pdf
  if (method_type=="pearson")
  {
    corr <- ccaPP::corPearson(as.vector(l_dmat), as.vector(dtree_dmat))
    
  } else if (method_type=="spearman") 
  {
    corr <- ccaPP::corSpearman(as.vector(l_dmat), as.vector(dtree_dmat), consistent = FALSE)
    
  } else if (method_type=="kendall")
  {
    corr <- ccaPP::corKendall(as.vector(l_dmat), as.vector(dtree_dmat), consistent = FALSE)

  } else
  {
    stop(paste0("Unsupported correlation method (mat_dtree_corr) requested"))
  }
  
}


matrix_item_size_check <- function(feature_length, dist_mat)
{
  ## check matrix size (they should have the same size)
  pre_matrix_size=0
  matrix_size_check_flag=T
  
  for (cy1 in seq(1:length(m_feature_length)))
  {
    matrix_size = length(as.vector(m_dist_mat[[cy1]]))
    
    if (pre_matrix_size!=0 && pre_matrix_size!=matrix_size)
    {
      print(paste0("l=", m_feature_length[[cy1]], "|", matrix_file_list[[cy1]], "|", matrix_size))
      matrix_size_check_flag=F
    }
    
    pre_matrix_size = matrix_size
  }
  
  if (matrix_size_check_flag==F) #stop going further if matrix sizes are vary
  {
    print("# input matrices sizes are vary, exit")
    stop("# input matrices sizes are vary")
    
  }  
}


#' @title visual_optimum_for_FFP
#' @Description Input distance matrices, calculate, organize, and output a dataframe for plotting purpose, in parallel computing environment
#' @param matrix_file_list path_to_distance_matrices (string)
#' @param tree_file_list path_to_phylogeny_tree (string)
#' @return dataframe for plotting FFP optimum
#' @examples  ffP_visual_optimum_parallel("./distance_matrices", "matrix"). Pick up only the files have 'matrix' string in names
optimum_determine_parallel <- function(m_combined_df, matrix_file_list, tree_file_list, tree_algorithm, outgroup_str, adjacent_flag, n_cores)
{
  #for sorting purpose; avoid lexigraphical order
  sfeature_length <- unlist(
    lapply(
      matrix_file_list
      , function(fn) as.numeric(tail(strsplit(fn, split="\\.")[[1]], n=1))
    )
  )
  
  matrix_file_list <- matrix_file_list[sort.list(sfeature_length)] # sort matrix file order in ascending order of feature length
  tree_file_list <- tree_file_list[sort.list(sfeature_length)] # sort matrix file order in ascending order of feature length
  sfeature_length <- sort(sfeature_length) #then sort sfeature_length
  
  m_combined_df["matrix_file_list"] <- matrix_file_list #error when dataframe is predefined of nrow=0
  
  ## run in parallel environment
  #https://julien-arino.github.io/2019/skel-parLapply
  # n_cores <- max(c(2, parallel::detectCores()-4)) #use four less cores, and minimum 2
  
  l_dist_mat <- parallel::mclapply( #load distance matrices
    X = matrix_file_list
    , function(x) read_dist_mat(x) #return symmetric distance matrix
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclapply
  )


  # l_dmat <- parallel::mclapply( #load distance matrices
  #   X = matrix_file_list
  #   , function(x) read_dist_mat(x) #return as.dist()
  #   , mc.cores = n_cores
  #   # , SIMPLIFY = F #unsupported in mclapply
  # )
  
  ## check matrix size (they should have the same size)
  matrix_item_size_check(sfeature_length, l_dist_mat)

  
  ## count the value >= upper_threshold
  dist_upper_threshold = 0.99999999 #equivalent to 1.0
  l_mat_thmax_count <- parallel::mclapply(
    X=l_dist_mat
    , fun = function(x) sum(as.vector(as.dist(x))>=dist_upper_threshold)
    , mc.cores = n_cores
  )
  
  
  l_emat_thmax_count <- list() #count values above 'dist_upper_threshold' among ingroups, after removing given outgroups
  if (outgroup_str!="none")
  {
    outgroup_list = unlist(strsplit(outgroup_str, ",")) #separated by commas
    
    l_emat_thmax_count <- parallel::mclapply(
      X=l_dist_mat
      , fun = function(x, y) sum(as.vector(as.dist(subset_dist_mat(x, y)))>=dist_upper_threshold)
      , y = outgroup_list
      , mc.cores = n_cores
    )
    
  } else
  {
    l_emat_thmax_count <- l_mat_thmax_count
  }

  
  # ret_list <- parallel::mclapply(
  #   X=load_path_list
  #   , function(x, y, z, f) hierarchical_cluster_func_annote(x, y, z, f)
  #   , y=annotation_table
  #   , z=save_folder_path
  #   , f=sqrt_flag
  #   , mc.cores=n_cores
  #   # , SIMPLIFY = F #unsupported in mclappy
  # )

  ## other matrix statistic
  l_dist_mat_median <- parallel::mclapply( #median
    X=l_dist_mat
    , fun=function(x) median(x)
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclappy
    )
  
  l_dist_mat_mean <- parallel::mclapply( #mean
    X=l_dist_mat
    , fun=function(x) mean(x)
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclappy
  )
  
  l_dist_mat_var <- parallel::mclapply( #varience(var) or standard deviation(sd)
    X=l_dist_mat
    , fun=function(x) var(x)
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclappy
  )
  
  ## now do multiprocessing using mcmapply (independent to cluster setup)
  # tree_algorithm="bionj"
  # corr_method = "kendall" #"pearson", "spearman", or "kendall"
  
  ## how to pass multiple arguments to parLapply (or lapply)  
  # https://stackoverflow.com/questions/14427253/passing-several-arguments-to-fun-of-lapply-and-others-apply
  # mclapply or lapply variants can accept X series plus addition arguments of fixed values (no serial like vector or list types)
  
  ## https://community.rstudio.com/t/use-multiple-arguments-in-mclapply/97404
  ## obtain distance-based phylo tree using specificed tree-building algorithm
  
  # subset dist() object https://rdrr.io/cran/usedist/man/dist_subset.html
  topo_shuffle_flag=T #shuffle matrix input order
  
  ## as.dist to phylotree
  l_phylo_tree <- parallel::mcmapply(function(l_dmat, tree_save_path, tree_algorithm, shuffle_flag, overwrite_flag)
    {build_phylo_tree(l_dmat, tree_save_path, tree_algorithm, shuffle_flag, overwrite_flag)}
    , l_dmat = l_dmat
    , tree_save_path = tree_file_list
    , tree_algorithm = rep(tree_algorithm, each=length(tree_file_list))
    , shuffle_flag= rep(shuffle_flag, each=length(tree_file_list))
    , overwrite_flag= rep(overwrite_flag, each=length(tree_file_list))
    , mc.cores = n_cores
    , SIMPLFIY = F
  )
  
  corr_method_types <- c("pearson", "spearman", "kendall") #three most popular correlation-r methods
  
  # calculate correlation-r between an input matrix and the tree-detrived matrix
  # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor
  ##calculate/obtain correlation-r of input_(JSD)matrix and the tree-derived distance matrix(input_matrix->phylo_tree->derived distance matrix)
  
  tree_dist_types <- c("rf_tree_dist") #c("rf_tree_dist", "spr_tree_dist") #do spr takes a lot of memory
  #skip spr_tree_dist fill with -1
  m_tree_dist_cm["spr_tree_dist"] = rep(-1, each=length(matrix_file_list))
  print(paste0("mcmapply skip, Adjacent ", "spr_tree_dist"))

  for (n_topology_distance in topo_dist_types)
  {
    m_combined_df[n_topology_distance] <- c(
      unlist(
        parallel::mcmapply(function(pri_l, lat_l, pri_tree, lat_tree, method_type, adjacent_flag) 
          {adjacent_topology_distance(pri_l, lat_l, pri_tree, lat_tree, method_type, adjacent_flag)}
          , pri_l = sfeature_length[1:(length(sfeature_length)-1)]
          , lat_l = sfeature_length[2:(length(sfeature_length))]
          , pri_tree = l_phylo_tree[1:(length(l_phylo_tree)-1)]
          , lat_tree = l_phylo_tree[2:(length(l_phylo_tree))]
          , method_type = rep(n_topology_distance, each=length(l_phylo_tree)-1)
          , adjacent_flag = rep(adjacent_flag, each=length(l_phylo_tree)-1)
          , mc.cores=n_cores
          , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
          )
        )
      , NA) # add NA at the end for no available comparison

    print(paste0("Calculated topology distance between two (adjacent) l-trees |", n_topology_distance))

  }

  ## Calculate correlation coefficient between two (adjacent) matrices
  for (n_corr_method in corr_method_types)
  {
    m_combined_df[paste0("mat_corr_", n_corr_method)] <- c(
      unlist(
        parallel::mcmapply(function(pri_l, lat_l, pri_dmat, lat_dmat, method_type, adjacent_flag)
          {adjacent_mat_corr(pri_l, lat_l, pri_dmat, lat_dmat, method_type, adjacent_flag)}
          , pri_l = sfeature_length[1:(length(sfeature_length)-1)]
          , lat_l = sfeature_length[2:(length(sfeature_length))]
          , pri_dmat = m_dist_mat[1:(length(l_dmat)-1)]
          , lat_dmat = m_dist_mat[2:(length(l_dmat))]
          , method_type = rep(n_corr_method, each=length(l_dmat)-1)
          , adjacent_flag = rep(adjacent_tight_flag, each=length(l_dmat)-1)
          , mc.cores=n_cores
          , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
      )
    )
    , NA) # add NA at the end for no available comparison
    
    print(paste0("Calculated correlation coefficient between two (adjacent) l-matrices |", n_corr_method))

  }

  ## Calculate correlation coefficient between input matrix and the tree-derived pairwise matrix
  ## this correlation coefficient indicates how well a phylotree represents the original distance matrix
  for (n_corr_method in corr_method_types)
  {
    m_combined_df[paste0("mat_corr_", n_corr_method)] <- unlist(
      parallel::mcmapply(function(l_dmat, l_tree, method_type)
        {mat_dtree_corr(l_dmat, l_tree, method_type)}
        , l_dmat = l_dmat
        , l_tree = l_phylo_tree
        , method_type = rep(n_corr_method, each=length(l_dmat))
        , mc.cores=n_cores
        , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
        )
    )
    
    print(paste0("Calculated tree representation correlation coefficient between a tree and the given matrix |", n_corr_method))
    
  }
  
  # stop()
  # combine all output into one dataframe
  m_combined_df <- data.frame(feature_length = unlist(m_feature_length)
                              , matrix_file_path = unlist(matrix_file_list)
                              , tree_file_path = unlist(tree_file_list)
                              
                              , matrix_median = unlist(m_matrix_median) #median
                              , matrix_mean = unlist(m_matrix_mean) #median
                              
                              , matrix_sd = unlist(m_matrix_sd) #standard deviation
                              , matrix_mad = unlist(m_matrix_mad) #median absolute deviation
                              
                              , matrix_upperlimit = unlist(m_matrix_upperlimit_count)
                              , matrix_wout_upperlimit = unlist(m_matrix_wout_upperlimit_count)
                              
                              # , phylo_tree = unlist(m_phylo_tree) #cannot coerce class 'phylo' in data.frame, unable
                              # , dtree_mat_corr_pearson = unlist(m_dtree_mat_corr_cm["dmcorr_pearson"]) #works but unintentionally assign row.names (instead of dmcorr_pearson)
                              , dtree_mat_corr_pearson = m_dtree_mat_corr_cm$dmcorr_pearson #unlist remove level to use designated column names (instead of dmcorr_pearson)
                              , dtree_mat_corr_spearman = m_dtree_mat_corr_cm$dmcorr_spearman
                              , dtree_mat_corr_kendall = m_dtree_mat_corr_cm$dmcorr_kendall
                              
                              , rf_tree_dist = m_tree_dist_cm$rf_tree_dist
                              , spr_tree_dist = m_tree_dist_cm$spr_tree_dist
                              
                              , mat_corr_pearson = m_mat_corr_cm$mcorr_pearson
                              , mat_corr_spearman = m_mat_corr_cm$mcorr_spearman
                              , mat_corr_kendall = m_mat_corr_cm$mcorr_kendall
                              
                              , stringsAsFactors = FALSE #as character
  )
  
  # row.names(m_combined_df) <- NULL #remove row.names when assigned unintentionally
  
  return (m_combined_df) #output to dataframe
}

#' @title manage_pre_data
#' @Description Read existing dataframe
#' @param file_path_list a list of files
#' @param df_load_path a working directory to load/save R dataframe
#' @param name_regex regular expression to specify files (string or NULL)
#' @return existing dataframe
#' @examples manage_pre_df(file_list, df_load_path, name_regex)
manage_pre_data <- function(current_wd, matrix_file_prefix=NULL, tree_file_prefix=NULL, outgroup_str, df_load_path) #list, path(str), bool
{
  
  # m_combined_df <- NULL
  # print(df_load_path)
  ## creat empty dataframe, instead of assigning NULL
  m_combined_df <- data.frame(feature_length = numeric()
                              , matrix_file_path = character()
                              , tree_file_path = character()
                              
                              , matrix_median = numeric()
                              , matrix_mean = numeric()
                              
                              , matrix_sd = numeric()
                              , matrix_mad = numeric()
                              
                              , matrix_upperlimit = numeric()
                              , matrix_wout_upperlimit = numeric()
                              
                              # , phylo_tree = unlist(m_phylo_tree) #cannot coerce class 'phylo' in data.frame, unable
                              , dtree_mat_corr_pearson = numeric()
                              , dtree_mat_corr_spearman = numeric()
                              , dtree_mat_corr_kendall = numeric()
                              
                              , rf_tree_dist = numeric()
                              , spr_tree_dist = numeric()
                              
                              , mat_corr_pearson = numeric()
                              , mat_corr_spearman = numeric()
                              , mat_corr_kendall = numeric()
                              
                              , stringsAsFactors = FALSE #as character
  )
  
  tree_algorithm = "bionj" #fixed tree algorithm
  
  template_colnames <- c()
  
  if (df_load_path=="none" || is.null(df_load_path))
  {
    
    return (m_combined_df) #NULL, nrow(df)==0
    
  } else if (file.exists(df_load_path))
  {
    
    print(paste0("read from: ", df_load_path))
    loaded_combined_df <- read.csv(df_load_path, header=T, comment.char="#")
    
    template_colnames <- colnames(m_combined_df)
    m_combined_df <- loaded_combined_df
    
    loaded_combined_df <- NULL #to free memory
    
  }
  
  print(paste0("current_wd: ", current_wd))
  setwd(current_wd) #change a working directory where matrix or tree files exist (advantageous when changing data folder); might cause continuous issue
  
  matrix_file_list <- noquote(
    list.files(
      path=current_wd
      , full.names=F #T: full file path, F: get only file name
      , pattern = paste0("^", matrix_file_prefix)
      , ignore.case=T
    )
  )
  
  ## to synchronize with version 4, specify a type of tree (nj, bionj, etc) in individual folders
  tree_folder = paste0(current_wd, "/", tree_algorithm)
  
  if (!file.exists(tree_folder))
  {
    dir.create(tree_folder)
    
  }
  
  #read tree files
  tree_file_list <- noquote(
    list.files(
      # path=current_wd
      path=tree_folder
      , full.names=F #T: full file path, F: get only file name
      , pattern = paste0("^", tree_file_prefix)
      , ignore.case=T
    )
  )
  
  # print(matrix_file_list)
  # print(tree_file_list)
  
  remove_m_feature <- c() #length 0
  add_m_feature <- c() #length 0
  
  print(paste("saved df ncol:", ncol(m_combined_df), "|", "required ncol(colnames):", length(template_colnames), sep=" "))
  
  # if (!is.null(m_combined_df) && ncol(m_combined_df)==df_item_cnt) #have saved df and meet the condition of a number of items
  if (nrow(m_combined_df)!=0 && (colnames(m_combined_df)==template_colnames)) #have saved df and meet the condition of a number of items
  {
    remove_m_feature <- setdiff(m_combined_df$matrix_file_path, matrix_file_list) #find new m_features not in saved df but in working directory
    add_m_feature <- setdiff(matrix_file_list, m_combined_df$matrix_file_path) #find new m_features not in saved df but in working directory
    
  } else #first time run (in case of m_combined_df==NULL)  # print(add_m_feature)
  {
    add_m_feature <- matrix_file_list
  }
  
  
  if (length(remove_m_feature)!=0) #since empty setdiff() = character(0) != NULL, but both lengths are 0
  {
    # print(paste0("remove: ", remove_m_feature))
    m_combined_df <- filter(m_combined_df, matrix_file_path!=remove_m_feature) #remove m_features only in saved df but not in working directory; not work when remove_m_feature==character(0)
    # m_combined_df <- subset(m_combined_df, matrix_file_path!=remove_m_feature) #remove m_features only in saved df but not in working directory
    # m_combined_df <- m_combined_df[m_combined_df$matrix_file_path!=remove_m_feature,] #subset row
    
  }
  
  
  if (length(add_m_feature)!=0 || (length(matrix_file_list)!=length(tree_file_list)))
  {
    print(paste0("add: ", add_m_feature))
    ## to synchronize with version 4, specify a type of tree (nj, bionj, etc) in individual folders
    
    if (length(matrix_file_list)!=length(tree_file_list))
    {
      tree_file_list <- paste0(tree_folder, "/", tree_file_prefix, tree_algorithm, ".", basename(matrix_file_list)) #assign tree file names    
    } else
    {
      tree_file_list <- paste0(tree_folder, "/", tree_file_list) #assign tree file names    
    }
    
    #print(paste0("tree:", tree_file_list))
    
    ## tree_file_list is dependent to matrix_file_list (update when new matrix files are added)
    # tree_file_list <- paste0(dirname(matrix_file_list), "/", tree_file_prefix, basename(matrix_file_list)) #use full filepath name
    # tree_file_list <- paste0(tree_file_prefix, basename(matrix_file_list)) #use file name
    
    m_combined_df <- ffp_visual_optimum_parallel(matrix_file_list, tree_file_list, outgroup_str, tree_algorithm)
    
  }
  
  
  if ((length(add_m_feature)!=0 || length(remove_m_feature)!=0) && !is.null(m_combined_df)) #save to a file only when there is any changes or contents
  {
    print(paste0("write to : ", df_load_path))
    write.csv(m_combined_df, df_load_path, row.names = FALSE, quote=FALSE) #omit writting to file for testing
    
  }
  
  # print(m_combined_df)
  return (m_combined_df)
}


opt <- parse_args(OptionParser(option_list=opt_list), positional_arguments = TRUE) #for trailing others

# standalone_run_all(project_list_load_path)

if (opt$option$all)
{
  print("check and run all")
  standalone_run_all(project_list_load_path)
  
} else if (opt$option$path!="")
{
  manage_pre_data(
    current_wd=opt$option$path
    , matrix_file_prefix=opt$option$matrixpre
    , tree_file_prefix=opt$option$treepre
    , outgroup_str=opt$option$outgroup
    , df_load_path=paste0(opt$option$path, ".df3")
  )
  
}


#' @title rf_contingency
#' @Description Find l-mer range of contingency
#' @param load_path path_to_distance_matrices (string)
#' @param name_regex regular expression to specify files (string or NULL)
#' @return dataframe for plotting FFP optimum
#' @examples  ffP_visual_optimum("./distance_matrices", "matrix"). Pick up only the files have 'matrix' string in names
rf_contingency <- function(lmer_array=c(), rf_array=c(), rf_min=0, n_elongation=2)
{
  #fine ranges of continuous RF minimum
  #using rle; output start and end l-mer ranges
  #https://stackoverflow.com/questions/43875716/find-start-and-end-positions-indices-of-runs-consecutive-values
  
  ##testing
  x <- c(1,2,3,5,6,7,8,0,0,0,0,1,12,3,4,5,0,0,0,0) #rf_array
  x.lmer <- c(2,3,4,5,6,7,8,9,10,15,16,17,18,30,31,32,50,51,52,53,54) #lmer_array
  
  rle_x <- rle(x)
  end_pos <- cumsum(rle_x$lengths)
  start_pos <- c(1, dplyr::lag(end_pos)[-1] + 1) #beware namespace conflict (lag)
    
  start_lmer <- x.lmer[start_pos]
  end_lmer <- x.lmer[end_pos]
  
  rf_contin <- data.frame(start=start_lmer, end=end_lmer, rf_value=rle_x$value, n_length=rle_x$lengths)
  # rf_contin
  
  rf_min=0
  n_elongation=3
  
  filter(rf_contin, rf_value == rf_min & n_length>=n_elongation) #considerate adjacent l-mer
  
  return(data.frame(start_pos, end_pos))
}


load_path = "/home/jjc/Desktop/ffp_optimum_plotly/Shiny_app/fungi_bionj_317/sym.mat.fp.20"
outgroup_list = c("rnd_min_1", "rnd_max_9")

load_path = "/home/jjc/Desktop/ffp_optimum_plotly/Shiny_app/primate_test_bionj/sym.mat.fp.24"
outgroup_list = c("9606_R")

ret_mat <- read_dist_mat(load_path)

upper_threshold=0.99999999

sum(as.vector(as.dist(subset_dist_mat(ret_mat, outgroup_list)))>=upper_threshold) #count excluding outgroups
sum(as.vector(as.dist(ret_mat))>=upper_threshold) #count including all

print(ret_mat)
ret_df <- as.data.frame(as.matrix(ret_mat))
colnames(ret_df)

