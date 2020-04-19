# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:02:16 2020
NOISYmputer python CORRECTING ERRORS
FILLING MISSING DATA
@author: Chrystian Sosa
@Correspondence author: Dr Mathias Lorieux
"""
#from numba import jit
import pandas as pd
import numpy as np
from numba import prange
import time
import copy

#https://www.pythonmania.net/es/2013/04/05/control-de-bucles-break-continue-y-pass/
#https://note.nkmk.me/en/python-numpy-where/
def chr_read_missing_data(chr_dir,n_chr):
       start_time = time.time()
       export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_5.txt")
       export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_6.txt")

       print("Starting step 6: %s" % export_file_path_filt)
       split_file = pd.read_csv(export_file_path_filt, sep=" ")
       chr_pos =  split_file["#"]
       sp_cols = list(split_file.columns)
       split_file2 = pd.DataFrame(columns = sp_cols)
       sp_cols.remove(sp_cols[0])
       split_file = split_file.drop(columns=['#'])
       split_file = split_file.to_numpy()
       print("Inds: %s" % str(split_file.shape[1]))
       print("Chr pos: %s" % str(split_file.shape[0]))
       
       for x in prange(split_file.shape[1]):
           #print("ind: %s" % str(x+1))
           tp=correct_loci_col_miss(split_file,x)
           split_file2.iloc[:,(x+1)] = tp

       split_file2["#"] = chr_pos
       split_file2.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
       split_file = pd.read_csv(export_file_path_filt, sep=" ")
       scl = changes_count(split_file2,split_file)
       print("--- %s seconds ---" % (time.time() - start_time))
       return(scl)

def correct_loci_col_miss(split_file,i):
    col_data = copy.deepcopy(split_file[:,i])
    col_length =  split_file.shape[0]
    j=0
    k=0
    allele_1 = ""
    allele_2 = ""
    while j < col_length-1:
        if col_data[j ]!="-":
            allele_1 = col_data[j]
        else:
            k = j+1
            while k < col_length:
                allele_2 = col_data[k]
                if allele_2 != "-":
                    if allele_1==allele_2:
                        sub2 = copy.deepcopy(col_data[j+1:k-1])
                        sub2=np.where(sub2=="-",allele_1 , sub2)
                        col_data[j+1:k-1]=sub2
                    else:
                        j=k+1
                        break
                k=k+1
                    #break
        j=j+1                         
    
    return(col_data)
    
    
    
def changes_count(split_file2,split_file):
    sp_cols = list(split_file2.columns)
    sp_cols.remove(sp_cols[0])
    counts_list=[]
    count_append = counts_list.append
    for a in range(len(sp_cols)):
        sp_count = split_file2[sp_cols[a]]==split_file[sp_cols[a]]
        sp_count = sp_count[sp_count==False].count()
        count_append(sp_count)
    s_cl = sum(counts_list)
    return(s_cl)    

chr_dir = "E:/CHR"  
#n_chr = 1
#chr_read_missing_data(chr_dir,n_chr)

x=[]
for a in prange(1,12+1):#12+1):
    n_chr = str(a)
    x.append(chr_read_missing_data(chr_dir,n_chr))