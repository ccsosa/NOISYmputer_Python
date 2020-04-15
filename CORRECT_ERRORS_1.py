# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:02:16 2020
NOISYmputer python CORRECTING ERRORS
CORRECT OBVIOUS ERRORS BEFORE FILLING 
WITH VERY HIGH THRESHOLD TO AVOID CORRECTING REAL HTZ CHUNKS
@author: Chrystian Sosa
@Correspondence author: Dr Mathias Lorieux
"""
#from numba import jit
import pandas as pd
import numpy as np
#from numba import jit, prange
import time
#from multiprocessing.pool import ThreadPool
#import subprocess
#@jit
#chr_dir,n_chr,popType

#@jit(parallel=True)
def chr_read(chr_dir,n_chr,popType,WindowSize):
    start_time = time.time()
    if popType == "SSD" or popType =="DH":
        export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_filtered_miss.txt")
        export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_filtered_step5.txt")
        print(export_file_path_filt+" Starting..."+ "| PopType="+popType+ " selected") 
       
        split_file = pd.read_csv(export_file_path_filt, sep=" ")
        chr_pos =  split_file["#"]
        sp_cols = list(split_file.columns)
        split_file2 = pd.DataFrame(columns = sp_cols)
        sp_cols.remove(sp_cols[0])
        split_file = split_file.drop(columns=['#'])
        split_file = split_file.to_numpy()
        print("Inds: "+ str(split_file.shape[1]) + " |" + " Chr pos:" + str(split_file.shape[0]))

        for i in range(split_file.shape[1]):
            #perc=np.round((i/((split_file.shape[1]))*100),3)
            #print(perc)
            tp=(correct_loci_col(split_file,WindowSize,sp_cols,i))
            split_file2.iloc[:,(i+1)] = tp
        
        split_file2["#"] = chr_pos
        split_file2.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
        split_file = pd.read_csv(export_file_path_filt, sep=" ")
        scl = changes_count(split_file2,split_file)
        print("--- %s seconds ---" % (time.time() - start_time))
        return(scl)
            
    elif popType == "F2" or popType == "BC1":
        print("--- %s popType ---" % (popType))
        print("--- %s seconds ---" % (time.time() - start_time))

        scl =0
        return(scl)

def correct_loci_col(split_file,WindowSize,sp_cols,i): 
    col_data = split_file[:,i]
    for j in range(len(split_file)):
        
        myRow1=j-WindowSize
        myRow2=j+WindowSize
        if myRow1 < 1:  myRow1 = 1
        if myRow2 > len(col_data): myRow2 = len(col_data)
        sub_col = col_data[[i for i in range(myRow1,myRow2) if i !=j]]
        sub_col = np.array_str(sub_col)

        nbA=sub_col.count("A")
        nbB=sub_col.count("B")
        nbH=sub_col.count("H")

        Nwindow = nbA + nbB + nbH
        if nbA >= Nwindow -1  and nbA > 4: #'we tolerate only one B or H
            if col_data[j]=="H":
                col_data[j] = "A"
        elif nbB >= Nwindow -1 and nbB > 4:
            if col_data[j]=="H":
                col_data[j] = "B"
        else:
            col_data[j] = col_data[j]
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
    #print("Sum of SNPs changed was", s_cl)
    
    
chr_dir = "E:/CHR"  
n_chr = str(1)
popType = "SSD"
WindowSize=7    
#num=6

x = chr_read(chr_dir,n_chr,popType,WindowSize)



    