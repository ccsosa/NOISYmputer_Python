# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 21:36:35 2020
IMPUTE SINGLETONS IN BREAKPOINTS
@author: cami_
"""
#def Nwindows_stats():
import pandas as pd
import numpy as np
from numba import prange
import time
import copy

def correct_singletons(chr_dir,n_chr,WindowSize): 
    start_time = time.time()
    export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_7.txt")
    export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_8.txt")
    print("Starting step 6: %s" % export_file_path_filt)
    split_file = pd.read_csv(export_file_path_filt, sep=" ")
    chr_pos =  split_file["#"]
    sp_cols = list(split_file.columns)
    split_file2 = pd.DataFrame(columns = sp_cols, index= split_file.index)
    sp_cols.remove(sp_cols[0])
    split_file = split_file.drop(columns=['#'])
    split_file = split_file.to_numpy()
    print("Inds: %s" % str(split_file.shape[1]))
    print("Chr pos: %s" % str(split_file.shape[0]))
 
    
    
def mark_single_rows(split_file,sp_cols,i):
        
    myAllele1=""
    myAllele2=""
    myAllele3=""
    st1= myAllele2=="A" and myAllele3=="B"  
    st2= myAllele2=="B" and myAllele3=="A"  
    st3= myAllele2=="B" and myAllele3=="H" 
    st4= myAllele2=="H" and myAllele3=="B"  
    st5= myAllele2=="A" and myAllele3=="H" 
    st6= myAllele2=="H" and myAllele3=="A"    
    
    col_data = split_file[:,i]
    col_data = copy.deepcopy(col_data)
    
    for j in range(1,len(col_data)-1):
        myAllele1= col_data[j]
        myAllele2= col_data[j-1]
        myAllele3= col_data[j+1]
        
        
        
        if('H' or '-') in myAllele1:
            if st1 or st2 or st3 or st4 or st5 or st6:
                if(j-(j-1)) < (j+1)-j:
                    mypos2 = mySNPPositions(rowSNP - 1)
                    myPos = mySNPPositions(rowSNP)
                    myPos3 = mySNPPositions(rowSNP + 1)
                    if myPos - mypos2 < myPos3 - myPos:
                        col_data[j]=myAllele2
                    elif myPos - mypos2 >= myPos3 - myPos:
                        col_data[j]=myAllele3
                        
            
                
    
    
    
    
    
    
    
chr_dir = "E:/CHR"  
n_chr = str(1)
WindowSize=30     
Chi2threshold = 3.84