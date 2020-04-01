# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:02:16 2020
NOISYmputer python FILTERING ON HTZ %, MISSING DATA %, MAF
@author: Chrystian Sosa
@Correspondence author: Dr Mathias Lorieux
"""
#IMPORTING MODULES
import pandas as pd
import numpy as np
import os.path
from numba import jit
    
nthreads = 4
    #size = 10**5  # CHANGED
@jit(parallel=True)
       
def filtering_snp_step(chr_dir,n_chr,maxFreqMD,maxFreqH,minFreqA,minFreqB):
    
        for a in range(1,n_chr+1):
            a1 = '{}'.format(a)
            #a1 = np.array2string(np.array(a))
            print("Chromosome "+ a1+" starting!")
            export_file_path = (chr_dir+"/"+"Chr"+a1+".txt")
            export_file_path_filt = (chr_dir+"/"+"Chr"+a1+"_filtered.txt")
            
            if(os.path.isfile(export_file_path_filt)):
                print("Chromosome "+ a1 + " already processed!")
                #split_file = type(None)
            else:
                split_file = pd.read_csv(export_file_path, sep=" ")
                #split_file = split_file.drop(columns=['chr'])
                nSNP = len(split_file)
                retain=[]
                
                for i in range(0,nSNP):
                    snp = split_file.iloc[i,:]
                    nbA = np.sum(snp.str.count('A'))
                    nbB = np.sum(snp.str.count('B'))
                    nbH = np.sum(snp.str.count('H'))
                    nbMD = np.sum(snp.str.count('-'))
                    Npop = nbA + nbB + nbH + nbMD
        
                    if(Npop > 0):
                        freqA = nbA / Npop
                        freqB = nbB / Npop
                        freqH = nbH / Npop
                        freqMD = nbMD / Npop
    
                        if(freqH < maxFreqH and freqMD < maxFreqMD and freqA > minFreqA and freqB > minFreqB):
                            #print("MATCH")
                            retain.append(True)
                        else:
                            #print("NOT MATCH")
                            retain.append(False)
                    else:
                        retain.append(False)
            
                split_file['filtered'] = retain
                split_file = split_file[split_file.filtered==1]
                split_file = split_file.drop(columns=['filtered'])
                split_file.to_csv(export_file_path_filt,index = False, header=True,sep=" ")
                print("Chromosome "+ a1+" done!")                
#        return(print("DONE"))
    
########################################3        


chr_dir = "E:/CHR"
maxFreqMD = 0.666
maxFreqH = 0.800
minFreqA = 0.010
minFreqB = 0.010
n_chr = 12
x_filtered = filtering_snp_step(chr_dir,n_chr,maxFreqMD,maxFreqH,minFreqA,minFreqB)
