# -*- coding: utf-8 -*-
"""
# -*- coding: utf-8 -*-
@author: cami_
Created on Thu Feb 27 11:02:16 2020
NOISYmputer python FILTERING ON HTZ %, MISSING DATA %, MAF
@author: Chrystian Sosa
@Correspondence author: Dr Mathias Lorieux
"""
#IMPORTING MODULES
import pandas as pd
import os.path
from tqdm import tqdm

#nthreads = 4
#size = 10**5  # CHANGED


def filtering_snp_step_chr(chr_dir,n_chrs,maxFreqMD,maxFreqH,minFreqA,minFreqB,log_dir):
    for i in tqdm(range(0,len(n_chrs)),desc="filtering snps"):
        n_chr=str(n_chrs[i])
        print(n_chr)
        filtering_snp_step(chr_dir,n_chr,maxFreqMD,maxFreqH,minFreqA,minFreqB)
    return("DONE!")

       
def filtering_snp_step(chr_dir,n_chr,maxFreqMD,maxFreqH,minFreqA,minFreqB):

    export_file_path = (chr_dir+"/"+"Chr"+n_chr+"_step_2.txt")
    export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_3.txt")

    retain =[]
    with open(export_file_path) as fp:
       for line in fp:
            nbA = line.count("A")
            nbB = line.count("B")
            nbH = line.count("H")
            nbMD = line.count("-")
            Npop = nbA + nbB + nbH + nbMD
            
            if Npop > 0:
                freqA = nbA / Npop
                freqB = nbB / Npop
                freqH = nbH / Npop
                freqMD = nbMD / Npop
    
                if freqH < maxFreqH and freqMD < maxFreqMD and freqA > minFreqA and freqB > minFreqB:
                    retain.append(True)
                else:
                    retain.append(False)
            else:
                retain.append(False)
    
    #retain.remove(0)
    del retain[0]
    if(os.path.isfile(export_file_path_filt)):
        return(print(export_file_path+" already processed!"))
    else:
        
        split_file2 = pd.read_csv(export_file_path, sep=" ")
        split_file = pd.read_csv(export_file_path, sep=" ")
        split_file['filtered'] = retain
        split_file = split_file[split_file.filtered==True]
        split_file = split_file.drop(columns=['filtered'])
#        split_file = split_file.drop_duplicates() 
        split_file.to_csv(export_file_path_filt,index = False, header=True,sep=" ")
         
       
       
        changed =split_file2.shape[0]-split_file.shape[0]
        data =[["step","prefiltering"],
                                  ["chr",n_chr],
                                  ["changes",changed]]
        log_file_df = pd.DataFrame(data,columns=["item","status"])
        log_file_df.to_csv(log_dir+"/"+"Chr"+n_chr+"_step3.log",index = False, header=True)
        
        return(print(export_file_path+" filtered!"))
########################################     


#chr_dir = "E:/CHR"
#maxFreqMD = 0.666
#maxFreqH = 0.800
#minFreqA = 0.010
#minFreqB = 0.010
#n_chrs=12
#n_chr = str(1)
#x_filtered = filtering_snp_step(chr_dir,n_chr,maxFreqMD,maxFreqH,minFreqA,minFreqB)

#filtering_snp_step_chr(chr_dir,n_chrs,maxFreqMD,maxFreqH,minFreqA,minFreqB)


