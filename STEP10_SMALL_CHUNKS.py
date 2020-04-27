# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:40:44 2020
DETECT IMPROBABLE SMALL CHUNKS
@author: cami_
"""
#from numba import jit
import pandas as pd
import numpy as np
from numba import prange
import time
import copy


#def small_chunks_step():    

   export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_9.txt")
   export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_10.txt")
       
   print("Starting: %s" % export_file_path_filt)
       
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
   


#def_SChunk_col():
i = 1
    col_data = split_file[:,i]
    col_data = copy.deepcopy(col_data)
    
    for j in range(1,len(col_data)))):
        myAllele = col_data[j]
        myAllele1 = col_data[j-1]
        if myAllele!=myAllele1:
            myStartRowChunk=j
            myStartPosChunk=chr_pos_num[j]
            k = j+1
            while k < len(col_data):
                if k> len(col_data):
                    break
                else:
                    myAllele2 =  col_data[k]
                    if myAllele2!=myAllele1:
                        myStopRowChunk = k
                        thisChunkSize = myStopRowChunk - myStartRowChunk + 1
                        chunkEnvironmentSize = (thisChunkSize * chunkEnvironmentMultiplier) + 1
                        if thisChunkSize <= smallChunkMaxSize:
                            myChunkCol = i
                            myStopPosChunk = chr_pos_num[myStopRowChunk]
                            chunkSizeBP = myStopPosChunk - myStartPosChunk
                            EnvironmentHomogeneous = True
                            startEnvirRow = myStartRowChunk - chunkEnvironmentSize
                            stopEnvirRow = myStopRowChunk + chunkEnvironmentSize
                            
                            if myAllele2 == myAllele1 or myAllele2 == "-" or myAllele1 == "-":
                                if myAllele1!="-": envirAllele = myAllele1
                                if myAllele2!="-": envirAllele = myAllele2
                                l = myStartRowChunk-1
                                while  l >= startEnvirRow-1:
                                    if l < 1: 
                                        EnvironmentHomogeneous = False
                                        break
                                    if  col_data[l]!=envirAllele:
                                        if  col_data[l]!="-":   
                                            EnvironmentHomogeneous = False
                                    l -= 1
                                    
                                l = myStopPosChunk+1   
                                while  l >= stopEnvirRow:
                                    if l >len(col_data): 
                                        EnvironmentHomogeneous = False
                                        break
                                    if  col_data[l]!=envirAllele:
                                        if  col_data[l]!="-":   
                                            EnvironmentHomogeneous = False
                                    l+= 1
                                    
                                if EnvironmentHomogeneous==True:
                                    startEnvirPos = chr_pos_num[startEnvirRow]
                                    stopEnvirPos = chr_pos_num[stopEnvirRow]
                                    EnvirSizeBP = stopEnvirPos - startEnvirPos
                                    if localMethod == 1:
                                        nbBKPTenvir = 0
                                        col_sub =np.array([split_file[j,m] for m in range(split_file.shape[1]) if m not in nbOfParents])                









smallChunkMaxSize = 100
chunkEnvironmentMultiplier = 0.55
localMethod = 2
minNbOfBKPTenvir = 6
minEnvirRecRate = 15.0
maxChunkRecRate = 50
nbOfParents= [0,1]# if it is 1 and 2 it should be [0,1] due to Python


smallChunkMaxSize 
chunkEnvironmentMultiplier
localMethod
minNbOfBKPTenvir
minEnvirRecRate
maxChunkRecRate
nbOfParents
chr_sep="_"
chr_dir = "E:/CHR"  
n_chr = str(1)
