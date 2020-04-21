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

def correct_singletons_break(chr_dir,n_chr,chr_sep): 
    start_time = time.time()
    export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_7.txt")
    export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_9.txt")
    print("Starting step 9: %s" % export_file_path_filt)
    split_file = pd.read_csv(export_file_path_filt, sep=" ")
    chr_pos =  split_file["#"]
    chr_pos_num = split_file["#"].str.split(chr_sep,n=1,expand=True)
    chr_pos_num = chr_pos_num.iloc[:,1]
    chr_pos_num = chr_pos_num.astype(int)
    sp_cols = list(split_file.columns)
    split_file2 = pd.DataFrame(columns = sp_cols, index= split_file.index)
    sp_cols.remove(sp_cols[0])
    split_file = split_file.drop(columns=['#'])
    split_file = split_file.to_numpy()
    print("Inds: %s" % str(split_file.shape[1]))
    print("Chr pos: %s" % str(split_file.shape[0]))
    for x in prange(split_file.shape[1]):
        #print("ind: %s" % str(x+1))
        tp=mark_single_rows(split_file,sp_cols,chr_pos_num,x)
        split_file2.iloc[:,(x+1)] = tp

    split_file2["#"] = chr_pos
    split_file2.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
    split_file = pd.read_csv(export_file_path_filt, sep=" ")
    scl = changes_count(split_file2,split_file)
    print("--- %s seconds ---" % (time.time() - start_time))
    return(scl)
    
def mark_single_rows(split_file,sp_cols,chr_pos_num,i):
    #i =1
    myAllele1=""
    myAllele2=""
    myAllele3=""
    st1= ""
    st2= ""  
    st3= "" 
    st4= ""  
    st5= ""
    st6= ""
    myPos=""
    mypos2=""
    myPos3=""
    
    col_data = split_file[:,i]
    col_data = copy.deepcopy(col_data)
    
    for j in range(1,len(col_data)-1):

        myAllele1= col_data[j]
        myAllele2= col_data[j-1]
        myAllele3= col_data[j+1]
        mypos2= chr_pos_num[j-1]
        myPos=chr_pos_num[j]
        myPos3=chr_pos_num[j+1]

        st1= myAllele2=="A" and myAllele3=="B"  
        st2= myAllele2=="B" and myAllele3=="A"  
        st3= myAllele2=="B" and myAllele3=="H" 
        st4= myAllele2=="H" and myAllele3=="B"  
        st5= myAllele2=="A" and myAllele3=="H" 
        st6= myAllele2=="H" and myAllele3=="A"
        
        if('H' or '-') in myAllele1:
            if st1 or st2 or st3 or st4 or st5 or st6:
                if myPos - mypos2 < myPos3 - myPos:
                    col_data[j] = myAllele2
                elif myPos - mypos2 >= myPos3 - myPos:
                    col_data[j] = myAllele3
        elif myAllele1 == "A":
            if st3 or st4:
                if myPos - mypos2 < myPos3 - myPos:
                    col_data[j] = myAllele2
                elif myPos - mypos2 >= myPos3 - myPos:
                    col_data[j] = myAllele3
                    
        elif myAllele1 == "B":
            if st5 or st6:  
                if myPos - mypos2 < myPos3 - myPos:
                    col_data[j] = myAllele2
                elif myPos - mypos2 >= myPos3 - myPos:
                    col_data[j] = myAllele3
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
n_chr = str(1)
chr_sep = "_"

x=[]
for a in prange(1,12+1):#12+1):
    n_chr = str(a)
    x.append(correct_singletons_break(chr_dir,n_chr,chr_sep))
