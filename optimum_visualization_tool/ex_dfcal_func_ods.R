# ---
# title: "Visually_determine_k-mer_optimum"
# author: "JaeJin Choi"
# date: "9/23/2020"
# ---

require(getopt)
require(dplyr)
require(stats)
require(ape) #contains phylo related functions; nj and bionj
require(phangorn) #contains upgma

# require(farver)
require(parallel) #for multi-core
# require(pcaPP) #https://rdrr.io/cran/pcaPP/man/cor.fk.html, for fast Kendall-tau estimation
require(ccaPP) #c++ implemented calculation; corKendall corPearson, corSpearman; https://search.r-project.org/CRAN/refmans/ccaPP/html/corFunctions.html
# require(cli) #load the latest version >=3.1.0

require(permute)

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



dist_matrix_sort <- function(input_mat) #necessary to input symmetric matrix (avoid triangular matrix that contains NAs)
{
  input_mat <- input_mat[order(rownames(input_mat)),]
  input_mat <- input_mat[, order(colnames(input_mat))]
  
  return(input_mat)
}



subset_dist <- function(dist_mat, exclude_list)
{
  sym_matrix <- as.matrix(dist_mat)
  
  ret_dist_matrix <- as.dist(sym_matrix[!(row.names(sym_matrix) %in% exclude_list) #row subset
                                        , !(colnames(sym_matrix) %in% exclude_list) #col subset
                                        ])
  
  return(ret_dist_matrix)
  
}



#' @title dtree_mat_corr
#' @Description Calculate correlation-r of input distance matrix and the tree-derived distance matrix
#' @param input_dist_mat distance_matrix
#' @param phylo_tree phylo tree object
#' @param corr_type specify correlation-r method (pearson, spearman, or kendall)
#' @return correlation-r
#' @examples read_to_dist_df(input_dist_mat, phylo_tree_obj, corr_type="")
dtree_mat_corr <- function(mdist, phylo_tree, corr_type="")
{
  ## https://rdrr.io/cran/ape/man/cophenetic.phylo.html; calculate tree nodes distances
  pw_tree_mdist <- as.dist(
    matrix_isort_sort(cophenetic.phylo(phylo_tree)) #get pairwise distance matrix from a phylogeny tree
    )
  
  corr <- 0.0
  ## http://www2.uaem.mx/r-mirror/web/packages/ccaPP/ccaPP.pdf
  if (corr_type=="pearson")
  {
    corr <- ccaPP::corPearson(as.vector(mdist), as.vector(pw_tree_mdist))
    
  } else if (corr_type=="spearman") 
  {
    corr <- ccaPP::corSpearman(as.vector(mdist), as.vector(pw_tree_mdist), consistent = FALSE)
    
  } else if (corr_type=="kendall")
  {
    corr <- ccaPP::corKendall(as.vector(mdist), as.vector(pw_tree_mdist), consistent = FALSE)

  }
  
  return(corr)

}


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


build_phylo_tree <- function(mdist, tree_save_path, tree_algorithm, shuffle_flag)
{
  output_tree <- ""
  
  if (shuffle_flag==TRUE)
  {
    rm_dist <- as.dist(matrix_iorder_shuffle(as.matrix(mdist)))
      
  }
  
  if (tree_algorithm=='bionj')
  {
    output_tree <- ape::bionj(dist_matrix_sh)
    
  } else if (tree_algorithm=='nj')
  {
    output_tree <- ape::nj(dist_matrix_sh)
    
  } else if (tree_algorithm=='upgma')
  {
    output_tree <- phangorn::upgma(dist_matrix_sh)
    
  }
  
  #write and save before output when does not exists
  write.tree(output_tree, file = tree_save_path)

  return (output_tree)
}


#' @title get_tree_dist
#' @Description Take a list of tree object and calculate their topology (Robinson-Foudls) differences in series 
#' @param tree_list a list of tree objects
#' @return A series of topology distance (between adjacent; n and n+1)
#' @examples  get_tree_dist(list_of_tree_objects)
get_tree_dist <- function(mphylo_trees)
{
  m_tree_dist <- as.numeric(c())
  #compare topology of adjacent trees (according to index)
  for (cy1 in seq(1:(length(mphylo_trees)-1))) #dist.topo(cy1, cy1+1)
  {
    # print(dist.topo(mphylo_trees[[cy1]], mphylo_trees[[cy1+1]]))
    m_tree_dist <- append(m_tree_dist, dist.topo(mphylo_trees[[cy1]], mphylo_trees[[cy1+1]])[1]) #just one element return; from list
    
  }
  
  m_tree_dist <- append(m_tree_dist, NA) #add last element as NA to fit size for m_combined_df
  
  return (m_tree_dist)
}


#' @title get_adjacent_input_compare_mcunit
#' @Description Take a list of tree object and calculate their topology (Robinson-Foudls) differences in series; using mcmapply 
#' @param pri_feature_length first l-mer
#' @param lat_feature_length second l-mer
#' @param pri_input first input
#' @param lat_input second input
#' @param comp_method_type specific method of comparison: Robinson-Foulds distance (tree_dist), and matrix correlation-r (pearson, spearman, or kendall)
#' @param adjacent_tight_flag only compare with adjacent l-mers (e.g., l-mer versus l-mer+1), else NA. Default is TRUE
#' @return output value
#' @examples  get_adjacent_input_compare_mcunit(pri_feature_length, lat_feature_length, pri_input, lat_input, comp_method_type="", adjacent_tight_flag=TRUE)
get_adjacent_input_compare_mcunit <- function(pri_feature_length, lat_feature_length, pri_input, lat_input, comp_method_type="", adjacent_tight_flag=TRUE) #except n_index, the rest given as constant paramenters
{
  output_value = NULL
  
  if (adjacent_tight_flag==TRUE && (lat_feature_length - pri_feature_length!=1))
  {
    output_value <- NA #is NA
    
  } else if (comp_method_type=="rf_tree_dist") #Robinson-Foulds(RF) topo-distance of two trees
  {
    #output_value <- dist.topo(pri_input, lat_input)[1] #from ape package
    output_value <- phangorn::RF.dist(pri_input, lat_input, check.labels=TRUE, normalize=F)
    
  } else if (comp_method_type=="spr_tree_dist") #Subtree Pruning and Re-grafting(SPR) distance of two trees
  {
    output_value <- phangorn::SPR.dist(pri_input, lat_input)
    
  } else if (comp_method_type=="pearson")
  {
    output_value <- ccaPP::corPearson(as.vector(pri_input), as.vector(lat_input))
    
  } else if (comp_method_type=="spearman")
  {
    output_value <- ccaPP::corSpearman(as.vector(pri_input), as.vector(lat_input), consistent = FALSE)
    
  } else if (comp_method_type=="kendall")
  {
    output_value <- ccaPP::corKendall(as.vector(pri_input), as.vector(lat_input), consistent = FALSE)
    
  }
  
  # using stat::cor which is slower implementation
  # } else if (comp_method_type=="pearson" || comp_method_type=="spearman") #calculate correlation-r of two matrices
  # {
  #   output_value <- cor(as.vector(pri_input), as.vector(lat_input), method=comp_method_type)
  #   # output_value <- cor(pri_input, lat_input, method=comp_method_type)
  #   
  # }  
  
  # output_vector <- append(output_vector, NA) #add last element as NA to fit size for m_combined_df
  
  return (output_value)
}

## adjacent topological distance
adjacent_compare_unit <- function(
  pri_feature_length, lat_feature_length #feature lengths
  , pri_dist_mat, lat_dist_mat #distance matrices
  , pri_phylo_tree, lat_phylo_tree #phylo trees
  , method_type #configure comparison type
  , adjacent_flag
  )
{
  output_value <- NULL
  
  if (adjacent_flag==TRUE && (lat_feature_length - pri_feature_length)!=1)
  {
    output_value <- NA #is NA
    
  } else if (method_type=="topo_rf") #Robinson-Fould distance
  {
    
  } else if (method_type=="topo_spr") #Subtree pruning and regrafting; much slower than calculating RF
  {
    
  } else if (method_type=="cbe_pearson")
  {
    
  } else if (method_type=="cbe_spearman")
  {
    
  } else if (method_type=="cbe_kendall")
  {
    
  }
  
}


optimum_determine_parallel <- function(matrix_file_list, tree_algorithm, outgroup_str, n_cores)
{
  #for sorting purpose; avoid lexigraphical order
  sfeature_length <- unlist(
    lapply(
      matrix_file_list
      , function(fn) as.numeric(tail(strsplit(fn, split="\\.")[[1]], n=1))
    )
  )
  
  matrix_file_list <- matrix_file_list[sort.list(m_feature_length)] # sort matrix file order in ascending order of feature length
  sfeature_length <- sort(sfeature_length) #then sort sfeature_length
  
  ## run in parallel environment
  #https://julien-arino.github.io/2019/skel-parLapply
  # n_cores <- max(c(2, parallel::detectCores()-4)) #use four less cores, and minimum 2
  n_cores <- 2 #fixed
  
  l_dist_mat <- parallel::mclapply( #load distance matrices
    X = matrix_file_list
    , function(x) read_dist_mat(x) #return symmetric distance matrix
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclapply
  )

  ## check matrix size (they should have the same size)
  pre_matrix_size=0
  matrix_size_check_flag=T
  
  for (cy1 in seq(1:length(m_feature_length)))
  {
    matrix_size = length(as.vector(m_dist_mat[[cy1]]))
    # print(paste0(cy1, "\t", matrix_size))
    
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
    # stop("# input matrices sizes are vary")
    return (NULL)
    
  }
  
  ## count the value >= upper_threshold
  l_mat_thmax_count <- parallel::mclapply(
    X=l_dist_mat
    , fun = function(x) sum(as.vector(as.dist(x))>=dist_upper_threshold)
    , mc.cores = n_cores
  )

  
  l_emat_thmax_count <- list()
  
  if (outgroup_str!="none")
  {
    outgroup_list = unlist(strsplit(outgroup_str, ",")) #separated by commas
    
    ## count the value >= upper_threshold, by excluding outgroups
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

  ## other matrix characteristics
  l_dist_mat_median <- parallel::mclapply(
    X=l_dist_mat
    , fun=function(x) median(x)
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclappy
    )
  
  l_dist_mat_mean <- parallel::mclapply(
    X=l_dist_mat
    , fun=function(x) mean(x)
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclappy
  )
  
  l_dist_mat_sd <- parallel::mclapply(
    X=l_dist_mat
    , fun=function(x) sd(x)
    , mc.cores = n_cores
    # , SIMPLIFY = F #unsupported in mclappy
  )
  
  ## exclude given artificial outgroups which are randomly generated through shuffling
  # subset dist() object https://rdrr.io/cran/usedist/man/dist_subset.html
  #print(dist_subset(m_dist_mat[[1]], c("10847", "T17"))) #works and can't drop but extract
  m_wout_dist_mat<-list() #empty list

  
  ## check size difference before and after exclusion
  # print(colnames(as.matrix(m_dist_mat[[1]])))
  # print(colnames(as.matrix(m_wout_dist_mat[[1]])))
  #quit()
  
  m_matrix_median <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) median(fn))
  print("parLapply done, median()")
  
  m_matrix_mean <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) mean(fn))
  print("parLapply done, mean()")
  
  ### calculate standard deviation of distance matrix
  m_matrix_sd <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) sd(fn))
  print("parLapply done, sd()")
  
  ### calculate median absolute deviation of distance matrix
  m_matrix_mad <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) mad(fn))
  print("parLapply done, mad()")
  
  
  ## obtain total count of values >= 0.999999 #0.999999, which is JSD upperlimit
  ### from a whole matrix, calculate JSD upperlimit
  mat_upper_threshold = 0.99999999
  m_matrix_upperlimit_count <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) sum(fn>=mat_upper_threshold)) #0.99999999 is almost equivalent to 1.0
  print(paste0("parLapply done, sum(); matrix_element >= ", mat_upper_threshold))
  
  ### after removing outgroups from a whole matrix, calculate JSD upperlimit(should be artifical outgroups that distances are 1.0)
  m_matrix_wout_upperlimit_count <- parallel::parLapply(cl=cl, X=m_wout_dist_mat, fun=function(fn) sum(fn>=mat_upper_threshold)) #0.99999999 is almost equivalent to 1.0
  print(paste0("parLapply done, sum(); matrix_element (without outgroups) >= ", mat_upper_threshold))

  
  ## now do multiprocessing using mcmapply (independent to cluster setup)
  # tree_algorithm="bionj"
  # corr_method = "kendall" #"pearson", "spearman", or "kendall"
  
  ## how to pass multiple arguments to parLapply (or lapply)  
  # https://stackoverflow.com/questions/14427253/passing-several-arguments-to-fun-of-lapply-and-others-apply
  # mclapply or lapply variants can accept X series plus addition arguments of fixed values (no serial like vector or list types)
  
  ## https://community.rstudio.com/t/use-multiple-arguments-in-mclapply/97404
  ## obtain distance-based phylo tree using specificed tree-building algorithm
  
  topo_shuffle_flag=T #shuffle matrix input order
  
  m_phylo_tree <- parallel::mcmapply(function(x, y, z, f) {get_phylogeny_tree(x, y, z, f)}
                                     , x=tree_file_list
                                     , y=m_dist_mat
                                     , z=rep(tree_algorithm, each=length(tree_file_list))
                                     , f=rep(topo_shuffle_flag, each=length(tree_file_list))
                                     , mc.cores=n_cores
                                     , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
  )
  
  print(paste0("mcmapply done, phylogeny_tree_output(): ", tree_algorithm))
  # print(m_phylo_tree)
  
  corr_method_types <- c("pearson", "spearman", "kendall") #three most popular correlation-r methods
  
  # calculate correlation-r between an input matrix and the tree-detrived matrix
  # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor
  ##calculate/obtain correlation-r of input_(JSD)matrix and the tree-derived distance matrix(input_matrix->phylo_tree->derived distance matrix)
  
  m_dtree_mat_corr_cm <- data.frame(matrix(nrow = length(m_dist_mat))) #empty dataframe
  
  for (n_corr_method in corr_method_types)
  {
    m_dtree_mat_corr <- parallel::mcmapply(function(x, y, z) {dtree_mat_corr(x, y, z)}
                                           , x=m_dist_mat
                                           , y=m_phylo_tree
                                           , z=rep(n_corr_method, each=length(m_dist_mat))
                                           , mc.cores=n_cores
                                           , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
    )
    
    m_dtree_mat_corr_cm[paste0("dmcorr_", n_corr_method)] <- unlist(m_dtree_mat_corr) #stack
    print(paste0("mcmapply done, dtree_mat_corr(): ", n_corr_method))
  }
  # print(m_dtree_mat_corr_cm)
  
  
  
  adjacent_tight_flag=TRUE #TRUE to calculate between l-mer and l-mer+1 only, else NA
  
  # calculate Robinson-Foulds(RF) distance between the trees of two adjacent l-mers
  # v_tree_dist <- get_adjacent_input_compare(m_phylo_tree, m_feature_length, comp_method_type="tree_dist", adjacent_tight_flag=adjacent_tight_flag)
  m_tree_dist_cm <- data.frame(matrix(nrow = length(m_dist_mat))) #empty dataframe
  
  tree_dist_types <- c("rf_tree_dist") #c("rf_tree_dist", "spr_tree_dist") #do spr takes a lot of memory
  #skip spr_tree_dist fill with -1
  m_tree_dist_cm["spr_tree_dist"] = rep(-1, each=length(matrix_file_list))
  print(paste0("mcmapply skip, Adjacent ", "spr_tree_dist"))
  
  #ape::dist.topo() and phangorn::RF.dist are equivalent = Robinson Foulds tree distance
  # v_spr_tree_dist <- get_adjacent_input_compare(m_phylo_tree, m_feature_length, comp_method_type="spr_tree_dist", adjacent_tight_flag=adjacent_tight_flag)
  # phangorn::sprdist() output an atomic vector of spr, spr_extra, rf, hdist (4 different outputs)
  # print(v_spr_tree_dist[["spr"]]) #get elements in an atomic vector using [[]] instead of $
  
  ## call specific tree distance only
  # print(phangorn::RF.dist(m_phylo_tree[[1]], m_phylo_tree[[2]], check.labels=TRUE))
  # print(phangorn::SPR.dist(m_phylo_tree[[1]], m_phylo_tree[[2]]))
  
  for (n_tree_dist_types in tree_dist_types)
  {
    m_tree_dist <- parallel::mcmapply(function(a, b, c, d, e, f) {get_adjacent_input_compare_mcunit(a, b, c, d, e, f)}
                                      , a = m_feature_length[1:(length(m_feature_length)-1)]
                                      , b = m_feature_length[2:(length(m_feature_length))]
                                      , c = m_phylo_tree[1:(length(m_phylo_tree)-1)]
                                      , d = m_phylo_tree[2:(length(m_phylo_tree))]
                                      
                                      , e = rep(n_tree_dist_types, each=length(m_phylo_tree)-1)
                                      , f = rep(adjacent_tight_flag, each=length(m_phylo_tree)-1)
                                      , mc.cores=n_cores
                                      , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
    )
    m_tree_dist <- c(m_tree_dist, NA) #for the last one; fill in with NA at the end
    print(paste0("mcmapply done, Adjacent ", n_tree_dist_types))
    
    m_tree_dist_cm[n_tree_dist_types] <- unlist(m_tree_dist) #stack
  }
  
  # calculate correlation-r between the matrices of two adjacent l-mers
  # v_mat_corr <- get_adjacent_input_compare(m_dist_mat, m_feature_length, comp_method_type=corr_method, adjacent_tight_flag=adjacent_tight_flag)
  # calculate correlation-r of three major methods
  m_mat_corr_cm <- data.frame(matrix(nrow = length(m_dist_mat))) #empty dataframe
  
  for (n_corr_method in corr_method_types)
  {
    m_mat_corr <- parallel::mcmapply(function(a, b, c, d, e, f) {get_adjacent_input_compare_mcunit(a, b, c, d, e, f)}
                                     , a = m_feature_length[1:(length(m_feature_length)-1)]
                                     , b = m_feature_length[2:(length(m_feature_length))]
                                     , c = m_dist_mat[1:(length(m_dist_mat)-1)]
                                     , d = m_dist_mat[2:(length(m_dist_mat))]
                                     , e = rep(n_corr_method, each=length(m_dist_mat)-1)
                                     , f = rep(adjacent_tight_flag, each=length(m_dist_mat)-1)
                                     , mc.cores=n_cores
                                     , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
    )
    
    m_mat_corr <- c(m_mat_corr, NA) #for the last one; fill in with NA at the end
    
    m_mat_corr_cm[paste0("mcorr_", n_corr_method)] <- unlist(m_mat_corr) #stack
    print(paste0("mcmapply done, Adjacent matrix corr done: ", n_corr_method))
  }
  # print(m_mat_corr_cm)
  
  
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




#' @title visual_optimum_for_FFP
#' @Description Input distance matrices, calculate, organize, and output a dataframe for plotting purpose, in parallel computing environment
#' @param matrix_file_list path_to_distance_matrices (string)
#' @param tree_file_list path_to_phylogeny_tree (string)
#' @return dataframe for plotting FFP optimum
#' @examples  ffP_visual_optimum_parallel("./distance_matrices", "matrix"). Pick up only the files have 'matrix' string in names
ffp_visual_optimum_parallel <- function(matrix_file_list, tree_file_list, outgroup_str, tree_algorithm) #current path, pattern=NULL
{
  #for sorting purpose; avoid lexigraphical order
  m_feature_length <- unlist(
    lapply(matrix_file_list
           , function(fn) as.numeric(tail(strsplit(fn, split="\\.")[[1]], n=1))
    )
  )
  
  # print(matrix_file_list)
  # print(tree_file_list)
  
  matrix_file_list <- matrix_file_list[sort.list(m_feature_length)] #sort file_list by ascending feature_length order
  tree_file_list <- tree_file_list[sort.list(m_feature_length)] #sort file_list by ascending feature_length order
  
  m_feature_length <- sort(m_feature_length) #also sort
  # print(m_feature_length)
  
  #https://julien-arino.github.io/2019/skel-parLapply
  # n_cores <- max(c(2, parallel::detectCores()-4)) #use four less cores, and minimum 2
  n_cores <- 2 #fixed
  
  #initiate cluster
  cl <- parallel::makeCluster(n_cores)
  print(paste0("parallel: initiate with #cores: ", n_cores))
  
  parallel::clusterExport(cl, c("read_matrix_to_df"
                                , "dist_matrix_sort" #sort by rownames and colnames, and then return as.dist()
                                # , "get_phylogeny_tree" #custom function function name plus two extra parameters
                                # , "get_adjacent_input_compare" #for adjacent l-mer comparison
                                , "sd", "mad", "median", "mean" #standard deviation, median absolute deviation, median; mean
                                , "bionj", "upgma", "nj" #tree building algorithms
                                , "read.tree", "write.tree" #ape; tree read and write
                                , "cophenetic.phylo" #calculate distances between terminal nodes
                                , "cor" #correlation coefficient
                                , "subset_dist"
  )
  , envir=environment()
  ) # Export needed variables (such as defined functions)
  
  print(paste0("parallel: export variables"))
  
  
  ## replace all lapply to parLapply that takes much time
  ## the rest lapply process in the element order of file_list
  m_dist_mat <- parallel::parLapply(cl=cl, X=matrix_file_list, fun=function(fn) read_matrix_to_df(fn)) #can accept either symmetric or lower-triangular matrix
  print("parLapply done, read_matrix_to_df()")
  
  ## check matrix size (they should have the same size)
  pre_matrix_size=0
  matrix_size_check_flag=T
  
  for (cy1 in seq(1:length(m_feature_length)))
  {
    matrix_size = length(as.vector(m_dist_mat[[cy1]]))
    # print(paste0(cy1, "\t", matrix_size))
    
    if (pre_matrix_size!=0 && pre_matrix_size!=matrix_size)
    {
      print(paste0("l=", m_feature_length[[cy1]], "|", matrix_file_list[[cy1]], "|", matrix_size))
      matrix_size_check_flag=F
    }
    
    pre_matrix_size = matrix_size
  }
  
  if (matrix_size_check_flag==F) #stop going further if matrix sizes are vary
  {
    parallel::stopCluster(cl)
    print("cluster object, release")
    print("# input matrices sizes are vary, exit")
    # stop("# input matrices sizes are vary")
    return (NULL)
    
  }
  
  
  ## exclude given artificial outgroups which are randomly generated through shuffling
  # subset dist() object https://rdrr.io/cran/usedist/man/dist_subset.html
  #print(dist_subset(m_dist_mat[[1]], c("10847", "T17"))) #works and can't drop but extract
  m_wout_dist_mat<-list() #empty list
  
  if (outgroup_str!="none")
  {
    outgroup_list = unlist(strsplit(outgroup_str, ",")) #separated by commas
    # print(outgroup_list)
    
    # m_wout_dist_mat <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) subset_dist(fun, exclude_list), exclude_list=outgroup_list) #not works
    m_wout_dist_mat <- parallel::parLapply(cl=cl, X=m_dist_mat, subset_dist, exclude_list=outgroup_list) #works
    
  } else
  {
    m_wout_dist_mat <- m_dist_mat
  }
  
  ## check size difference before and after exclusion
  # print(colnames(as.matrix(m_dist_mat[[1]])))
  # print(colnames(as.matrix(m_wout_dist_mat[[1]])))
  #quit()
  
  m_matrix_median <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) median(fn))
  print("parLapply done, median()")
  
  m_matrix_mean <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) mean(fn))
  print("parLapply done, mean()")
  
  ### calculate standard deviation of distance matrix
  m_matrix_sd <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) sd(fn))
  print("parLapply done, sd()")
  
  ### calculate median absolute deviation of distance matrix
  m_matrix_mad <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) mad(fn))
  print("parLapply done, mad()")
  
  
  ## obtain total count of values >= 0.999999 #0.999999, which is JSD upperlimit
  ### from a whole matrix, calculate JSD upperlimit
  mat_upper_threshold = 0.99999999
  m_matrix_upperlimit_count <- parallel::parLapply(cl=cl, X=m_dist_mat, fun=function(fn) sum(fn>=mat_upper_threshold)) #0.99999999 is almost equivalent to 1.0
  print(paste0("parLapply done, sum(); matrix_element >= ", mat_upper_threshold))
  
  ### after removing outgroups from a whole matrix, calculate JSD upperlimit(should be artifical outgroups that distances are 1.0)
  m_matrix_wout_upperlimit_count <- parallel::parLapply(cl=cl, X=m_wout_dist_mat, fun=function(fn) sum(fn>=mat_upper_threshold)) #0.99999999 is almost equivalent to 1.0
  print(paste0("parLapply done, sum(); matrix_element (without outgroups) >= ", mat_upper_threshold))
  
  
  # Stop cluster. no more parallel compute
  parallel::stopCluster(cl)
  print("cluster object, release")
  
  ## now do multiprocessing using mcmapply (independent to cluster setup)
  # tree_algorithm="bionj"
  # corr_method = "kendall" #"pearson", "spearman", or "kendall"
  
  ## how to pass multiple arguments to parLapply (or lapply)  
  # https://stackoverflow.com/questions/14427253/passing-several-arguments-to-fun-of-lapply-and-others-apply
  # mclapply or lapply variants can accept X series plus addition arguments of fixed values (no serial like vector or list types)
  
  ## https://community.rstudio.com/t/use-multiple-arguments-in-mclapply/97404
  ## obtain distance-based phylo tree using specificed tree-building algorithm
  
  topo_shuffle_flag=T #shuffle matrix input order
  
  m_phylo_tree <- parallel::mcmapply(function(x, y, z, f) {get_phylogeny_tree(x, y, z, f)}
                                     , x=tree_file_list
                                     , y=m_dist_mat
                                     , z=rep(tree_algorithm, each=length(tree_file_list))
                                     , f=rep(topo_shuffle_flag, each=length(tree_file_list))
                                     , mc.cores=n_cores
                                     , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
  )
  
  print(paste0("mcmapply done, phylogeny_tree_output(): ", tree_algorithm))
  # print(m_phylo_tree)
  
  corr_method_types <- c("pearson", "spearman", "kendall") #three most popular correlation-r methods
  
  # calculate correlation-r between an input matrix and the tree-detrived matrix
  # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor
  ##calculate/obtain correlation-r of input_(JSD)matrix and the tree-derived distance matrix(input_matrix->phylo_tree->derived distance matrix)
  
  m_dtree_mat_corr_cm <- data.frame(matrix(nrow = length(m_dist_mat))) #empty dataframe
  
  for (n_corr_method in corr_method_types)
  {
    m_dtree_mat_corr <- parallel::mcmapply(function(x, y, z) {dtree_mat_corr(x, y, z)}
                                           , x=m_dist_mat
                                           , y=m_phylo_tree
                                           , z=rep(n_corr_method, each=length(m_dist_mat))
                                           , mc.cores=n_cores
                                           , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
    )
    
    m_dtree_mat_corr_cm[paste0("dmcorr_", n_corr_method)] <- unlist(m_dtree_mat_corr) #stack
    print(paste0("mcmapply done, dtree_mat_corr(): ", n_corr_method))
  }
  # print(m_dtree_mat_corr_cm)
  
  
  
  adjacent_tight_flag=TRUE #TRUE to calculate between l-mer and l-mer+1 only, else NA
  
  # calculate Robinson-Foulds(RF) distance between the trees of two adjacent l-mers
  # v_tree_dist <- get_adjacent_input_compare(m_phylo_tree, m_feature_length, comp_method_type="tree_dist", adjacent_tight_flag=adjacent_tight_flag)
  m_tree_dist_cm <- data.frame(matrix(nrow = length(m_dist_mat))) #empty dataframe
  
  tree_dist_types <- c("rf_tree_dist") #c("rf_tree_dist", "spr_tree_dist") #do spr takes a lot of memory
  #skip spr_tree_dist fill with -1
  m_tree_dist_cm["spr_tree_dist"] = rep(-1, each=length(matrix_file_list))
  print(paste0("mcmapply skip, Adjacent ", "spr_tree_dist"))
  
  #ape::dist.topo() and phangorn::RF.dist are equivalent = Robinson Foulds tree distance
  # v_spr_tree_dist <- get_adjacent_input_compare(m_phylo_tree, m_feature_length, comp_method_type="spr_tree_dist", adjacent_tight_flag=adjacent_tight_flag)
  # phangorn::sprdist() output an atomic vector of spr, spr_extra, rf, hdist (4 different outputs)
  # print(v_spr_tree_dist[["spr"]]) #get elements in an atomic vector using [[]] instead of $
  
  ## call specific tree distance only
  # print(phangorn::RF.dist(m_phylo_tree[[1]], m_phylo_tree[[2]], check.labels=TRUE))
  # print(phangorn::SPR.dist(m_phylo_tree[[1]], m_phylo_tree[[2]]))
  
  for (n_tree_dist_types in tree_dist_types)
  {
    m_tree_dist <- parallel::mcmapply(function(a, b, c, d, e, f) {get_adjacent_input_compare_mcunit(a, b, c, d, e, f)}
                                      , a = m_feature_length[1:(length(m_feature_length)-1)]
                                      , b = m_feature_length[2:(length(m_feature_length))]
                                      , c = m_phylo_tree[1:(length(m_phylo_tree)-1)]
                                      , d = m_phylo_tree[2:(length(m_phylo_tree))]
                                      
                                      , e = rep(n_tree_dist_types, each=length(m_phylo_tree)-1)
                                      , f = rep(adjacent_tight_flag, each=length(m_phylo_tree)-1)
                                      , mc.cores=n_cores
                                      , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
    )
    m_tree_dist <- c(m_tree_dist, NA) #for the last one; fill in with NA at the end
    print(paste0("mcmapply done, Adjacent ", n_tree_dist_types))
    
    m_tree_dist_cm[n_tree_dist_types] <- unlist(m_tree_dist) #stack
  }
  
  # calculate correlation-r between the matrices of two adjacent l-mers
  # v_mat_corr <- get_adjacent_input_compare(m_dist_mat, m_feature_length, comp_method_type=corr_method, adjacent_tight_flag=adjacent_tight_flag)
  # calculate correlation-r of three major methods
  m_mat_corr_cm <- data.frame(matrix(nrow = length(m_dist_mat))) #empty dataframe
  
  for (n_corr_method in corr_method_types)
  {
    m_mat_corr <- parallel::mcmapply(function(a, b, c, d, e, f) {get_adjacent_input_compare_mcunit(a, b, c, d, e, f)}
                                     , a = m_feature_length[1:(length(m_feature_length)-1)]
                                     , b = m_feature_length[2:(length(m_feature_length))]
                                     , c = m_dist_mat[1:(length(m_dist_mat)-1)]
                                     , d = m_dist_mat[2:(length(m_dist_mat))]
                                     , e = rep(n_corr_method, each=length(m_dist_mat)-1)
                                     , f = rep(adjacent_tight_flag, each=length(m_dist_mat)-1)
                                     , mc.cores=n_cores
                                     , SIMPLIFY = F #to return as a list, SIMPLIFY = T as a default to reduce dimension
    )
    
    m_mat_corr <- c(m_mat_corr, NA) #for the last one; fill in with NA at the end
    
    m_mat_corr_cm[paste0("mcorr_", n_corr_method)] <- unlist(m_mat_corr) #stack
    print(paste0("mcmapply done, Adjacent matrix corr done: ", n_corr_method))
  }
  # print(m_mat_corr_cm)
  
  
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

