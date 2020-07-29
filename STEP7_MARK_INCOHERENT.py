# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 21:47:37 2020
 MARK AND FILTER INCOHERENT LOCI
@author: cami_
"""

    

#def Nwindows_stats():
import pandas as pd
import numpy as np
import copy
from tqdm import tqdm

def correct_loci_row_chr(chr_dir,n_chrs,Chi2threshold,WindowSize,log_dir): 
    for i in tqdm(range(len(n_chrs)),desc="Mark and filter incoherent loci"):
        n_chr=str(n_chrs[i])
        print(n_chr)
        correct_loci_row(chr_dir,n_chr,Chi2threshold,WindowSize,log_dir)
        
    return("DONE!")
    
def correct_loci_row(chr_dir,n_chr,Chi2threshold,WindowSize,log_dir): 
    export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_6.txt")
    export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_7.txt")
    split_file = pd.read_csv(export_file_path_filt, sep=" ")
    chr_pos =  split_file["#"]
    sp_cols = list(split_file.columns)
    split_file2 = pd.DataFrame(columns = sp_cols, index= split_file.index)
    sp_cols.remove(sp_cols[0])
    split_file = split_file.drop(columns=['#'])
    split_file = split_file.to_numpy()
    
    chi2SNP1=[]
    chi2SNP2=[]
    NbOfBadSNP=[]
    
#    ExpA=0;ExpH=0;ExpB=0;Npop_j=0;Npop=0
    for j in range(split_file.shape[0]):
        row_data = split_file[j,:]
        row_data = copy.deepcopy(row_data)
        
        myRow1=j-WindowSize
        myRow2=j+WindowSize
        if myRow1 < 1:  myRow1 = 1
        if myRow2 > split_file.shape[0]: myRow2 = split_file.shape[0]
        sub_col = split_file[[i for i in range(myRow1,myRow2) if i !=j],:]
        nbA=float(len(sub_col[sub_col=="A"]))
        nbB=float(len(sub_col[sub_col=="B"]))
        nbH=float(len(sub_col[sub_col=="H"]))
        #nbMD=float(sub_col.count("-"))
        Npop = float(nbA + nbB + nbH)
        ExpA = nbA / Npop; ExpB = nbB / Npop;# ExpH = nbH / Npop
###########################################
        sub_col_j = copy.deepcopy(split_file[j,:])
        nbA_j=float(len(sub_col_j[sub_col_j=="A"]))
        nbB_j=float(len(sub_col_j[sub_col_j=="B"]))
        nbH_j=float(len(sub_col_j[sub_col_j=="H"]))
          #nbMD_j=float(sub_col_j.count("-"))
        Npop_j = nbA_j + nbB_j + nbH_j
        
################        
        if (nbA_j - ExpA * Npop_j)**2 ==0 or (ExpA * Npop_j) ==0:
            exp1=0
        else:
            exp1=(nbA_j - ExpA * Npop_j)  **2 / (ExpA * Npop_j) 
            
################
        if (nbB_j - ExpB * Npop_j)** 2  ==0 or (ExpB * Npop_j) ==0:
            exp2=0
        else:
            exp2=(nbB_j - ExpB * Npop_j)** 2 / (ExpB * Npop_j)
            
        x2N1=exp1+exp2
        chi2SNP1.append(x2N1)
################################################################
################################################################
################################################################        
        if (nbB_j - ExpA * Npop_j) ** 2  ==0 or (ExpA * Npop_j) ==0:
            exp1a=0
        else:
            exp1a=(nbB_j - ExpA * Npop_j) ** 2 / (ExpA * Npop_j)
            ################
        if (nbA_j - ExpB * Npop_j)  **2 ==0 or (ExpB * Npop_j) ==0:
            exp2a=0
        else:
            exp2a=(nbA_j - ExpB * Npop_j)  **2 / (ExpB * Npop_j) 

        x2N2=exp1a+exp2a
        chi2SNP2.append(x2N2)    
################################################################
        if x2N1>Chi2threshold:
            if x2N2>Chi2threshold:
                 
                NbOfBadSNP.append(j)
                chi2SNP2.append(1)
                row_data=row_data
            else:
                #NbOfBadSNP.append(str(j)+":fix")
                sub2 = copy.deepcopy(row_data)
                sub2=np.where(sub2=="A","X" , sub2)
                sub2=np.where(sub2=="B","Y" , sub2)
                sub2=np.where(sub2=="X","B" , sub2)
                sub2=np.where(sub2=="Y","A" , sub2)
                row_data = sub2
        else: 
            row_data=row_data
        split_file2.iloc[j,range(1,split_file2.shape[1])]=  row_data

    split_file = pd.read_csv(export_file_path_filt, sep=" ")
    if(len(NbOfBadSNP)>0):
        print(f'removing {len(NbOfBadSNP)} loci!')
        split_file2 = split_file2.drop(NbOfBadSNP)
    else:
        print(f'removing {0} loci!')
        
    split_file2["#"] = chr_pos
    split_file2.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")


    changes=(split_file.shape[0]-split_file2.shape[0])*split_file2.shape[1]
    data =[["step","Mark and filter incoherent loci"],
           ["chr",n_chr],
           ["NbOfBadSNP",len(NbOfBadSNP)],
           ["changes",changes]]
    
    log_file_df = pd.DataFrame(data,columns=["item","status"])
    log_file_df.to_csv(log_dir+"/"+"Chr"+n_chr+"_step5.log",index = False, header=True)