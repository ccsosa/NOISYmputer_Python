# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:34:20 2020

@author: cami_
"""

#IMPORTING MODULES
import sys
import pandas as pd
from tqdm import tqdm

def remove_dup_snp_step_chr(chr_dir,n_chrs,missing_opt):
    
    for i in tqdm(range(1,n_chrs+1),desc="filtering snps"):
        n_chr=str(i)
        remove_dup_snp_step(chr_dir,n_chr,missing_opt=missing_opt)
    return("DONE!")
    

def remove_dup_snp_step(chr_dir,n_chr,missing_opt):

    export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_3.txt")
    export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_4.txt")
    split_file = pd.read_csv(export_file_path_filt, sep=" ")
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
        return(print(export_file_path_filt2+" filtered!"))
    elif(missing_opt=="classic"):
        print("using missing data")
        split_file = split_file.drop_duplicates(subset=sp_cols)
        split_file.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
        return(print(export_file_path_filt2+" filtered!"))
    else:
        print("error,not valid method chosen")
        sys.exit()
    
    
   


chr_dir = "E:/CHR"
#n_chr = str(1)

remove_dup_snp_step_chr(chr_dir,n_chrs,missing_opt="missing")