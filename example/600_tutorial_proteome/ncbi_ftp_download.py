# NCBI FTP download
# FTP address: ftp.ncbi.nlm.nih.gov
# FTP FAQ: http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/

import os, sys, getopt, shutil
from ftplib import FTP

import gzip, io

from Bio import SeqIO

from multiprocessing import cpu_count, Pool

import pandas as pd

import time


EXCLUDE_FASTA_KEYWORDS_LIST=[ #keywords to exclude; organelle related keywords
    '(mitochondrion)'
    , '(plasmid)'
    , '(chloroplast)'
    , 'mitochondrion, complete genome'
    , 'mitochondrial'
    , 'chloroplast, complete genome'
    , 'plast, complete genome'
    ]


def multi_process(proc_name, args_list=[], n_cpu=0): #args_list is a list of tuple arguments set

    ret_list=[]

    with Pool(processes=n_cpu) as proc_pool:
        ret_proc = proc_pool.starmap_async(proc_name, args_list, chunksize=int(len(args_list) / n_cpu))

        proc_pool.close()
        proc_pool.join()

        ret_list = ret_proc.get()


    return(ret_list)


def ftp_download_unit(
    ftp_address=""
    , ftp_path=""
    , save_path=""
    , exclude_fasta_keywords_list=[]
    , download_flag=False
    ):

    total_fasta_item_cnt=0
    filtered_fasta_item_cnt=0
    
    file_name = os.path.basename(ftp_path)

    ftp = FTP(ftp_address)
    ftp.login()

    ftp.cwd(os.path.dirname(ftp_path)) #change path, necessary to use ftp.nlst()

    bfile_io = io.BytesIO()
    str_io = io.StringIO()


    if file_name in ftp.nlst():
        ftp.retrbinary(f"RETR {file_name}", bfile_io.write) ##download file
    
    ftp.quit()

    
    if (bfile_io.tell()==0): #file exists but empty
        print(f"# Either not found or is empty: {ftp_address}/{ftp_path}")

    else:
        bfile_io.seek(0)
        
        try:
            with gzip.open(bfile_io, "rt") as gread_f: #"rt" -> read as text
                for seq_ob in SeqIO.parse(gread_f, "fasta"):
        
                    total_fasta_item_cnt+=1
                    fasta_header = f">{seq_ob.description}"
                    
                    if (True not in [True if seq_ob.description.lower().find(n_keyword)!=-1 else False for n_keyword in exclude_fasta_keywords_list]):
                        
                        if (download_flag==True):
                            str_io.write(f"{fasta_header}\n") #str_io.write(f">{seq_ob.description}\n")
                            str_io.write(f"{str(seq_ob.seq)}\n")
                            str_io.write("\n") #add line break between items for read friendly
                            
                        filtered_fasta_item_cnt+=1
        
        except Exception as err:
            print(repr(err))
    
        if (str_io.tell()!=0 and (os.path.exists(save_path)==False or str_io.tell()!=os.path.getsize(save_path))):
            str_io.seek(0)
            
            with open(save_path, 'w') as write_f:
                shutil.copyfileobj(str_io, write_f)
                
    str_io.truncate(0)
    str_io.close()
    
    bfile_io.truncate(0)
    bfile_io.close()
    

    return([f"{ftp_address}/{ftp_path}", total_fasta_item_cnt, filtered_fasta_item_cnt])



def show_help():
    print('Usage: python [program_path] [options] [file1] [file2] ...')
    print('Options:')
    print('-h, show_help')
    print('-S [path], save_folder_path, create if not exists. Current working directory is default')
    print('-p [str], primary_column_name. Default is "taxid"')
    print('-a [str], attribute_column_name. Default is "NCBI FTP download links"')
    print('-t [int], # of thread to run at once. Default is # of CPU - 1')
    print('-d, do download. Default is off')
    print('-e [str], --exclude_keywords, comma separated list of keywords to exclude. Default is organelle related keywords')
    print(f"\t{', '.join(EXCLUDE_FASTA_KEYWORDS_LIST)}")

    print("\n")
    print("# After the run it will summarize and print a log, either download flag is on or off`")

    sys.exit()


if __name__=='__main__':

    opts, args = getopt.getopt(sys.argv[1:], 'hS:p:a:t:de:', ['thread=', "exclude_keywords="])
    
    save_folder_path = os.getcwd() #='/home/jjc/Desktop' #for test
    ftp_address = 'ftp.ncbi.nlm.nih.gov' #default

    primary_colname = "taxid"
    attribute_colname = "NCBI FTP download links"
    
    download_flag=False

    n_cpu = cpu_count() - 1

    exclude_fasta_keywords_list=EXCLUDE_FASTA_KEYWORDS_LIST
  
    for opt, arg in opts:
        
        if (opt=='-h'):
            show_help()

        elif (opt=='-S'):
            save_folder_path=os.path.abspath(arg)

        elif (opt=='-p'):
            primary_colname=str(arg)

        elif (opt=='-a'):
            attribute_colname=str(arg)

        elif (opt=='-t' or opt=='--thread'): # a number of threads
            n_cpu=max(1, int(arg)) #least one
        
        elif (opt=='-d'):
            download_flag=True

        elif (opt=='-e' and opt=='--exclude_keywords'):
            exclude_fasta_keywords_list = [n_item.strip().lower() for n_item in arg.split(',')]


    combine_df=pd.DataFrame()

    for n_arg in args:
        load_path = os.path.abspath(n_arg)

        ''' #for tab separated text file
        pd_df = pd.read_csv(
            filepath_or_buffer=load_path
            , comment="#" #skip the line begin with comment (single char)
            , header=0
            , skip_blank_lines=True
            , sep="\t"
            , na_values=['na', 'NA', 'none', 'None', 'N//A', 'n//a', '', 'nan', '-']
            , dtype=str
            )
        '''

        pd_df = pd.read_excel(
            io=load_path
            , comment="#" #skip the line begin with comment (single char)
            , header=0
            #, skip_blank_lines=True
            #, sep="\t"
            , na_values=['na', 'NA', 'none', 'None', 'N//A', 'n//a', '', 'nan', '-']
            , dtype=str
            )

        if (combine_df.empty):
            combine_df = pd_df

        else:
            combine_df = pd.concat(
                [combine_df, pd_df]
                , ignore_index=True
                , join="inner"
                , verify_integrity=True #check if the index is unique
                )

    args_list=[]

    for n_name, n_ftp_file in combine_df[[primary_colname, attribute_colname]].itertuples(index=False):

        ftp_address = n_ftp_file.split('/')[0]
        ftp_path = n_ftp_file.replace(ftp_address, "")
        save_path = f"{save_folder_path}/{n_name}"

        args_list.append(
            [ftp_address, ftp_path, save_path, exclude_fasta_keywords_list, download_flag]
        )


    if (len(args_list)>0):

        if (os.path.exists(save_folder_path)==False):
            os.makedirs(save_folder_path)

        ret_list = multi_process(ftp_download_unit, args_list, n_cpu)

        pd_df = pd.DataFrame(
            ret_list
            , columns=[attribute_colname, "total_fasta_item_cnt", "filtered_fasta_item_cnt"]
            )
        #print(pd_df.columns)
        #print(combine_df.columns)

        combine_df = combine_df.join(
            pd_df.set_index(attribute_colname)
            , on=attribute_colname
            , how="left"
            , rsuffix="_dup"
            )
        
        combine_df.drop(
            [i for i in combine_df.columns if '_dup' in i]
            , axis=1
            , inplace=True
        )

        print(f"# Run result, obtained / submitted = {len(pd_df[pd_df['filtered_fasta_item_cnt'] > 0])} / {len(pd_df)}")
        print(f"## excluded/preprocessed by excluding FASTA header (and the sequence) that contains keywords,")
        print(f"##\t{'; '.join(exclude_fasta_keywords_list)}")
        print(f"## success if 'filtered_fasta_item_cnt' > 0")

        log_summary_path = f"{save_folder_path}.report.log"
        
        with open(log_summary_path, "w") as write_f:
            write_f.write(f"# Summarized at {time.strftime("%Y-%m-%d", time.gmtime())}\n") #date, import time
            write_f.write(f"# Run result, obtained / submitted = {len(pd_df[pd_df['filtered_fasta_item_cnt'] > 0])} / {len(pd_df)}\n")
            write_f.write(f"## excluded/preprocessed by excluding FASTA header (and the sequence) that contains keywords,\n")
            write_f.write(f"##\t{'; '.join(exclude_fasta_keywords_list)}\n")
            write_f.write(f"## success if 'filtered_fasta_item_cnt' > 0\n")
            
            combine_df.to_csv(
                path_or_buf=write_f #log_summary_path
                , sep="\t"
                , index=False
                , header=True
            )
            
        print(f'# Log saved at {log_summary_path}')
        print("")

    print("# Work done")
