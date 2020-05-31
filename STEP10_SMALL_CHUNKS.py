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


#def small_chunks_step(me):    

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
        myAllele1 = col_data[j-1] #allele just before the candidate chunk
        if myAllele!=myAllele1: #possible start of a chunk detected
            myStartRowChunk=j #store start position
            myStartPosChunk=chr_pos_num[j] #in BP
            k = j+1
            while k < len(col_data):
                if k> len(col_data):
                    break
                else:
                    myAllele2 =  col_data[k] #allele of SNP just after end of candidate chunk
                    if myAllele2!=myAllele1: #possible end of chunk detected
                        myStopRowChunk = k
                        thisChunkSize = myStopRowChunk - myStartRowChunk + 1 # chunk size in number of rows
                        chunkEnvironmentSize = (thisChunkSize * chunkEnvironmentMultiplier) + 1
                        if thisChunkSize <= smallChunkMaxSize: #it's a chunk of acceptable size
                            myChunkCol = i
                            myStopPosChunk = chr_pos_num[myStopRowChunk] #in BP
                            chunkSizeBP = myStopPosChunk - myStartPosChunk
                            EnvironmentHomogeneous = True
                            startEnvirRow = myStartRowChunk - chunkEnvironmentSize
                            stopEnvirRow = myStopRowChunk + chunkEnvironmentSize
                            
                            if myAllele2 == myAllele1 or myAllele2 == "-" or myAllele1 == "-": #verify that the environment is homogeneous
                                if myAllele1!="-": envirAllele = myAllele1
                                if myAllele2!="-": envirAllele = myAllele2
                                l = myStartRowChunk-1 
                                while  l >= startEnvirRow-1:#verify homogeneity up
                                    if l < 1: 
                                        EnvironmentHomogeneous = False
                                        break
                                    if  col_data[l]!=envirAllele:
                                        if  col_data[l]!="-":   #some MD can be around
                                            EnvironmentHomogeneous = False
                                    l -= 1
                                    
                                l = myStopPosChunk+1   
                                while  l >= stopEnvirRow: #verify homogeneity down
                                    if l >len(col_data): 
                                        EnvironmentHomogeneous = False
                                        break
                                    if  col_data[l]!=envirAllele: #not homogeneous, stop
                                        if  col_data[l]!="-":   
                                            EnvironmentHomogeneous = False
                                    l+= 1
                                    
                                if EnvironmentHomogeneous==True:
                                    #'now we need the local probability of recombination to decide if we replace alleles
                                    startEnvirPos = chr_pos_num[startEnvirRow]
                                    stopEnvirPos = chr_pos_num[stopEnvirRow]
                                    EnvirSizeBP = stopEnvirPos - startEnvirPos #envir size in bp

                                    EnvircMSize = 0




if mappingFunction == "Kosambi": 
    'conversion to centimorgans
    F1cM = Kosambi(F1RecFreq)
elif mappingFunction == "Haldane":
    F1cM = Haldane(F1RecFreq)
elif mappingFunction = "none" :
    F1cM = 100 * F1RecFreq

def Kosambi(Recomb_Fraction):

    if Recomb_Fraction == -1:
        Kosambi=-1
    elif Recomb_Fraction == 0:
        Kosambi = 0
    elif Recomb_Fraction < 0.499:
        Kosambi = 25 * np.log((1 + 2 * Recomb_Fraction) / (1 - 2 * Recomb_Fraction))
    elif Recomb_Fraction >= 0.499:
        Kosambi = 25 * np.log((1 + 2 * 0.499) / (1 - 2 * 0.499))
    return(Kosambi)

def Haldane(Recomb_Fraction):

    if Recomb_Fraction == -1:
        Haldane=-1
    elif Recomb_Fraction == 0:
        Haldane = 0
    elif Recomb_Fraction < 0.499:
        Haldane = -50 * np.log(1 - 2 * Recomb_Fraction)
    elif Recomb_Fraction >= 0.499:
        Haldane = -50 * np.log(1 - 2 * 0.499)
    return(Haldane)






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


"""
                                    if localMethod == 1:
                                        
                                        this method shouldn't be used
                                        'method 1, just count recombinants in the chunk environment
                                        'count number of breakpoints in environment
                                        
                                        nbBKPTenvir = 0
                                        col_sub =np.array([split_file[j,m] for m in range(split_file.shape[1]) if m not in nbOfParents])                
                                        for x in range(startEnvirRow,stopEnvirRow):
                                            if col_sub[x]!=col_sub[x-1]:
                                                nbBKPTenvir += 1
                                         
                                            """