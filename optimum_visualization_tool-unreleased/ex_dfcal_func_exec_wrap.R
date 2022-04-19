# ---
# title: "Visually_determine_k-mer_optimum"
# author: "JaeJin Choi"
# date: "9/23/2020"
# ---

## load "dfcal_func.R" and execute without web-visualization

require(optparse) #R optget in python-style
opt_list = list(
  make_option(c("-a", "--all"), type="character", action="store_true", default=F, help="Check and run all projects enlisted in project_load_list (.ods or .txt)")
  , make_option(c("-p", "--path"), type="character", action="store", default="", help="Set specific project working directory")
  , make_option(c("-m", "--matrixpre"), type="character", action="store", default="mat.", help="Set matrix file pattern")
  , make_option(c("-t", "--treepre"), type="character", action="store", default="tree.", help="Set tree file pattern")
  , make_option(c("-o", "--outgroup"), type="character", action="store", default="none", help="Set outgroup(s), delimit with comma','")
  # , make_option(c("-p", "--path"), type="character", action="store", default="", help="Set specific project working directory")
  
  # make_option(c("-m", "--maintitle"), type="character", action="store", default="", help="Set plot's main title"),
  # make_option(c("-s", "--subtitle"), type="character", action="store", default="", help="Set plot's subtitle"),
  # make_option(c("-c", "--caption"), type="character", action="store", default="", help="Set plot's caption/description underneath figure"),
  # make_option(c("-S", "--savefile_path"), type="character", action="store", default="", help="Set savefile name"),
  # make_option(c("-k", "--kit"), type="integer", action="store", default=0, help="k=1 (RF plot; smoothing), k=2 (distanec distribution), k=3 (line graph)"),
  # make_option(c("-e", "--opttemplate"), type="character", action="store", default="", help="value template")
)


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
