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
from numba import prange
import time
#import multiprocessing as mp
import copy
#from multiprocessing import Pool
#from multiprocessing.pool import ThreadPool
#import subprocess
#@jit
#chr_dir,n_chr,popType
#https://numba.pydata.org/numba-doc/latest/user/parallel.html
#https://medium.com/@mjschillawski/quick-and-easy-parallelization-in-python-32cb9027e490

def chr_read(chr_dir,n_chr,popType,WindowSize):
    start_time = time.time()
    if popType == "SSD" or popType =="DH":
        export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_4.txt")
        export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_5.txt")
        
        print("Starting: %s" % export_file_path_filt)
        print("PopType: %s" % popType)
       
        split_file = pd.read_csv(export_file_path_filt, sep=" ")
        chr_pos =  split_file["#"]
        sp_cols = list(split_file.columns)
        split_file2 = pd.DataFrame(columns = sp_cols)
        sp_cols.remove(sp_cols[0])
        split_file = split_file.drop(columns=['#'])
        split_file = split_file.to_numpy()
        print("Inds: %s" % str(split_file.shape[1]))
        print("Chr pos: %s" % str(split_file.shape[0]))
        
        start_time = time.time()
        
        for i in prange(split_file.shape[1]):
            tp=correct_loci_col(split_file,WindowSize,sp_cols,i)
            split_file2.iloc[:,(i+1)] = tp
         
        split_file2["#"] = chr_pos
        split_file2.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
        split_file = pd.read_csv(export_file_path_filt, sep=" ")
        scl = changes_count(split_file2,split_file)
        print("--- %s seconds ---" % (time.time() - start_time))
        return(scl)
            
    elif popType == "F2" or popType == "BC1":
       print("PopType: %s" % popType)
       scl =0
       return(scl)
        

def correct_loci_col(split_file,WindowSize,sp_cols,i): 
    col_data = split_file[:,i]
    col_data = copy.deepcopy(col_data)
    for j in range(len(col_data)):
        
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
#n_chr = str(1)
popType = "SSD"
WindowSize=7    
#num=mp.cpu_count()
#print("Number of processors: ", mp.cpu_count())
#x = chr_read(chr_dir,n_chr,popType,WindowSize,num)

x=[]
for a in prange(1,12+1):
    nchr = str(a)
    x.append(chr_read(chr_dir,nchr,popType,WindowSize))

"""
x=[]
for a in prange(1,12+1):
    nchr = str(a)
    x.append(chr_read(chr_dir,nchr,popType,WindowSize,num))

pool = mp.Pool(mp.cpu_count())
results = pool.starmap(chr_read, [(chr_dir,n_chr,popType,WindowSize,num) for n_chr in list(range(1,13))])
pool.close()
tp2 = ThreadPool(num)
for a in range(1,12+1):
    nchr = str(a)        
    tp2.apply_async(chr_read(chr_dir,nchr,popType,WindowSize,num))
tp2.close()
"""


