# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:34:20 2020

@author: cami_
"""

#IMPORTING MODULES
import sys
import pandas as pd
from tqdm import tqdm

def remove_dup_snp_step_chr(chr_dir,n_chrs,missing_opt,log_dir):
    
    for i in tqdm(range(len(n_chrs)),desc="filtering snps"):
        n_chr=str(n_chrs[i])
        print(n_chr)
        remove_dup_snp_step(chr_dir,n_chr,missing_opt,log_dir)
    return("DONE!")
    

def remove_dup_snp_step(chr_dir,n_chr,missing_opt,log_dir):

    export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_3.txt")
    export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_4.txt")
    split_file = pd.read_csv(export_file_path_filt, sep=" ")
    split_file2 = pd.read_csv(export_file_path_filt, sep=" ")

    sp_cols = list(split_file.columns)
    sp_cols.remove(sp_cols[0])
    
    if(missing_opt=="missing"):
        print("skipping missing data")
        ListToStr_list= []
        for i in range(0,len(split_file)):
            listToStr = [str(x) for x in split_file.iloc[i,1:]]
            listToStr = ' '.join(map(str, listToStr))
            listToStr = listToStr.translate({ord(i): None for i in '-'})
            listToStr = listToStr.translate({ord(' '): None})
            ListToStr_list.append(listToStr)
#          ListToStr_list.remove(0)
        split_file['STRING'] = ListToStr_list
        split_file = split_file.drop_duplicates('STRING')
        split_file = split_file.drop(columns=['STRING'])
        split_file.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
    elif(missing_opt=="classic"):
        print("using missing data")
        split_file = split_file.drop_duplicates(subset=sp_cols)
        split_file.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
    else:
        print("error,not valid method chosen")
        sys.exit()
    
    
    changed =split_file2.shape[0]-split_file.shape[0]
    data =[["step","prefiltering"],
                                  ["chr",n_chr],
                                  ["changes",changed]]
    log_file_df = pd.DataFrame(data,columns=["item","status"])
    log_file_df.to_csv(log_dir+"/"+"Chr"+n_chr+"_step4.log",index = False, header=True)
        
    return(print(export_file_path_filt2+" filtered!"))


#remove_dup_snp_step_chr(chr_dir,n_chrs,missing_opt="missing")