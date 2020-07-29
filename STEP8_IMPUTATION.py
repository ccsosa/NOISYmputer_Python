# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 21:36:35 2020
IMPUTATIOn
@author: cami_
"""
import pandas as pd
import numpy as np
from numba import prange
import copy
from tqdm import tqdm

def imputation_chr(chr_dir,chr_sep,n_chrs,chr_dir2,imputation_params_options,imputation_params_step1_params,imputation_params_step2_params,imputation_params_step3_params,log_dir): 
    
    
    for i in tqdm(range(len(n_chrs)),desc="Breakpoint singletons"):
        n_chr=str(n_chrs[i])
        imputation(chr_dir,chr_sep,n_chr,chr_dir2,imputation_params_options,imputation_params_step1_params,imputation_params_step2_params,imputation_params_step3_params,log_dir)
        
    return("DONE!")


def imputation(chr_dir,chr_sep,n_chr,chr_dir2,imputation_params_options,imputation_params_step1_params,imputation_params_step2_params,imputation_params_step3_params,log_dir): 
    export_file_path_filt = (chr_dir+"/"+"Chr"+n_chr+"_step_7.txt")
    export_file_path_filt2 = (chr_dir2+"/"+"Chr"+n_chr+"_step_8.txt")
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
    
    
    
    n_changes = []
    for x in prange(split_file.shape[1]):
        #print("ind: %s" % str(x+1))
        tp=imputation_col(split_file,imputation_params_options,imputation_params_step1_params,imputation_params_step2_params,imputation_params_step3_params,x)
        split_file2.iloc[:,(x+1)] = tp[0]
        n_changes.append(tp[1])
    split_file2["#"] = chr_pos
    split_file2.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
    split_file = pd.read_csv(export_file_path_filt, sep=" ")
 
    data =[["step","Imputation"],
           ["chr",n_chr],
           ["changes",sum(n_changes)]]
    log_file_df = pd.DataFrame(data,columns=["item","status"])
    log_file_df.to_csv(log_dir+"/"+"Chr"+n_chr+"_step8.log",index = False, header=True)


def imputation_col(split_file,imputation_params_options,imputation_params_step1_params,imputation_params_step2_params,imputation_params_step3_params,i):
    n_changes = []
    col_data = split_file[:,i]
    col_data = copy.deepcopy(col_data)
    x0 = imputation_0(split_file,i)
    n_changes.append(x0[1])
    
  
     #PARAMS METHOD 1
    WindowSize1=imputation_params_step1_params[0][1]
    minFreqA1=imputation_params_step1_params[1][1]
    minFreqB1=imputation_params_step1_params[2][1]
    minFreqBinH1=imputation_params_step1_params[3][1]
    minFreqAinH1=imputation_params_step1_params[4][1]
    HTZtrans_rate1=imputation_params_step1_params[5][1]
    
    #PARAMS METHOD 2
    WindowSize2=imputation_params_step2_params[0][1]
    minFreqA2=imputation_params_step2_params[1][1]
    minFreqB2=imputation_params_step2_params[2][1]
    minFreqBinH2=imputation_params_step2_params[3][1]
    minFreqAinH2=imputation_params_step2_params[4][1]
    
    #PARAMS METHOD 3
    WindowSize3=imputation_params_step3_params[0][1]

    #METHOD 1
    if imputation_params_options[0][1]==True & imputation_params_options[1][1]==False & imputation_params_options[2][1]==False:
        x_imp = imputation_met1(x0[0],HTZtrans_rate1,minFreqBinH1,minFreqAinH1,minFreqA1,minFreqB1,WindowSize1)
        n_changes.append(x_imp[1])
    #METHOD 2
    elif imputation_params_options[0][1]==False & imputation_params_options[1][1]==True & imputation_params_options[2][1]==False:
        x_imp = imputation_met2(x0[0],minFreqBinH2,minFreqAinH2,minFreqA2,minFreqB2,WindowSize2)
        n_changes.append(x_imp[1])
    #METHOD 3
    elif imputation_params_options[0][1]==False & imputation_params_options[1][1]==False & imputation_params_options[2][1]==True:
        x_imp = imputation_met3(x0[0],WindowSize3)
        n_changes.append(x_imp[1])
    #METHOD 1 &  METHOD 2
    elif imputation_params_options[0][1]==True & imputation_params_options[1][1]==True & imputation_params_options[2][1]==False:
        x_imp = imputation_met1(x0[0],HTZtrans_rate1,minFreqBinH1,minFreqAinH1,minFreqA1,minFreqB1,WindowSize1)
        n_changes.append(x_imp[1])
        x_imp = imputation_met2(x_imp[0],minFreqBinH2,minFreqAinH2,minFreqA2,minFreqB2,WindowSize2)
        n_changes.append(x_imp[1])
    #METHOD 1 &  METHOD 3
    elif imputation_params_options[0][1]==True & imputation_params_options[1][1]==False & imputation_params_options[2][1]==True:
        x_imp = imputation_met1(x0[0],HTZtrans_rate1,minFreqBinH1,minFreqAinH1,minFreqA1,minFreqB1,WindowSize1)
        n_changes.append(x_imp[1])
        x_imp = imputation_met3(x_imp[0],WindowSize3)
        n_changes.append(x_imp[1])
    #METHOD 2 & METHOD3
    elif imputation_params_options[0][1]==False & imputation_params_options[1][1]==True & imputation_params_options[2][1]==True:
         x_imp = imputation_met2(x0[0],minFreqBinH2,minFreqAinH2,minFreqA2,minFreqB2,WindowSize2)
         n_changes.append(x_imp[1])
         x_imp = imputation_met3(x_imp[0],WindowSize3)
         n_changes.append(x_imp[1])
    #METHOD 1, 2, and 3
    elif imputation_params_options[0][1]==True & imputation_params_options[1][1]==True & imputation_params_options[2][1]==True:
         x_imp = imputation_met1(x0[0],HTZtrans_rate1,minFreqBinH1,minFreqAinH1,minFreqA1,minFreqB1,WindowSize1)
         n_changes.append(x_imp[1])
         x_imp = imputation_met2(x_imp[0],minFreqBinH2,minFreqAinH2,minFreqA2,minFreqB2,WindowSize2)
         n_changes.append(x_imp[1])
         x_imp = imputation_met3(x_imp[0],WindowSize3)
         n_changes.append(x_imp[1])
         
    return(col_data,sum(n_changes))

def imputation_0(split_file,i):
    col_data = split_file[:,i]
    col_data = copy.deepcopy(col_data)
    NbSNP_Imputed =0
    for myRow in range(2,len(col_data)-2):
        #'compare SNP in the middle of the window with the window major allele  
        rowSNP = myRow
        myAllele = col_data[rowSNP]
        myAllele1 = col_data[rowSNP-2]
        myAllele2 = col_data[rowSNP-1]
        myAllele3 = col_data[rowSNP+1]
        myAllele4 = col_data[rowSNP+2]
        
        if myAllele != myAllele1:
            if myAllele1==myAllele2 and myAllele2 ==myAllele3 and myAllele3 == myAllele4:
                col_data[rowSNP] = myAllele1
                NbSNP_Imputed += 1 
    return(col_data,NbSNP_Imputed)
                
def imputation_met1(col_data,HTZtrans_rate,minFreqBinH,minFreqAinH,minFreqA,minFreqB,WindowSize):
      HTZtrans = int(WindowSize * 2 + 1) * HTZtrans_rate
      NbSNP_Imputed =0
      for myRow in range(2,len(col_data)-2):
        #'compare SNP in the middle of the window with the window major allele  
        rowSNP = myRow
        myAllele = col_data[rowSNP]
        myRow1 = rowSNP - WindowSize #'start of window
        myRow2 = rowSNP + WindowSize  #'end of window
        
        if myRow1 < 1:
            myRow1 = 1
            myRow2 = rowSNP + (rowSNP - 1) #'to avoid unbalanced windows (otherwise it suppresses small terminal chunks)
            if myRow2 > len(col_data):
                myRow2 = len(col_data)
                myRow1 = rowSNP - (len(col_data) - rowSNP)
                #'count transitions in window
                nbTrans = 0

                AlleleMemory = col_data[myRow1] #'read first allele
                for thisRow  in range(myRow1,myRow2):
                    if col_data[thisRow] != AlleleMemory:  #'might be a transition
                        if col_data[thisRow] != "-" and AlleleMemory != "-":# Then 'it's a true transition
                            nbTrans += 1
                            AlleleMemory = col_data[thisRow]
                            
                            
                #'allele frequencies
                nbA = 0; nbB = 0; nbH = 0; nbMD = 0; Nwindow = 0
                NbSNP_Imputed=0
                for thisRow  in range(myRow1,myRow2):
                    if col_data[thisRow] == "A":
                        nbA = nbA + 1
                    elif col_data[thisRow] == "B":
                        nbB = nbB + 1
                    elif col_data[thisRow] == "H":
                        nbH = nbH + 1
                    elif col_data[thisRow] == "-":
                        nbMD = nbMD + 1
                Nwindow = nbA + nbB + nbH + nbMD
                
                if Nwindow == 0: #this can happen when many consecutive missing data
                    col_data[rowSNP] = myAllele
                else:
                    freqA = nbA / Nwindow; freqB = nbB / Nwindow;# freqH = nbH / Nwindow
                    #'compare SNP in the middle of the window with the window major allele
                    col_data[rowSNP] = myAllele
                    
                    if freqA >= minFreqA:#'it should be an A
                        if myAllele != "A":
                            col_data[rowSNP] = "A"
                            NbSNP_Imputed += 1
                    elif freqB >= minFreqB:# Then 'it should be a B
                        if myAllele != "B":
                            col_data[rowSNP] = "B"
                            NbSNP_Imputed += 1
                    elif nbTrans >= HTZtrans:
                        if freqB >= minFreqBinH and freqA >= minFreqAinH:
                            if myAllele != "H":
                                col_data[rowSNP] = "H"
                            NbSNP_Imputed += 1

        return(col_data,NbSNP_Imputed)

def imputation_met2(col_data,maxFreqBinH,minFreqAinH,minFreqA,minFreqB,WindowSize):
    #problem it suppresses small terminal chunks  
    NbSNP_Imputed = 0
    for myRow in range(2,len(col_data)-2):
        #'compare SNP in the middle of the window with the window major allele  
        rowSNP = myRow
        myAllele = col_data[rowSNP]
        myRow1 = rowSNP - WindowSize #'start of window
        myRow2 = rowSNP + WindowSize  #'end of window
        
        if myRow1 < 1:
            myRow1 = 1
            myRow2 = rowSNP + (rowSNP - 1) #'to avoid unbalanced windows (otherwise it suppresses small terminal chunks)
            if myRow2 > len(col_data):
                myRow2 = len(col_data)
                myRow1 = rowSNP - (len(col_data) - rowSNP)
                #'calculate frequencies in the windownbA = 0; nbB = 0; nbH = 0; nbMD = 0; Nwindow = 0
                nbA = 0; nbB = 0; nbH = 0; nbMD = 0; Nwindow = 0
                

                for thisRow  in range(myRow1,myRow2):
                    if col_data[thisRow] == "A":
                        nbA = nbA + 1
                    elif col_data[thisRow] == "B":
                        nbB = nbB + 1
                    elif col_data[thisRow] == "H":
                        nbH = nbH + 1
                    elif col_data[thisRow] == "-":
                        nbMD = nbMD + 1
                Nwindow = nbA + nbB + nbH + nbMD
                if Nwindow == 0: #this can happen when many consecutive missing data
                    col_data[rowSNP] = myAllele
                else:
                    freqA = nbA / Nwindow; freqB = nbB / Nwindow; freqH = nbH / Nwindow
                    #'compare SNP in the middle of the window with the window major allele
                    col_data[rowSNP] = myAllele
                    
                    if freqA >= minFreqA:#'it should be an A
                        if myAllele != "A":
                            col_data[rowSNP] = "A"
                            NbSNP_Imputed += 1
                    elif freqB >= minFreqB:# Then 'it should be a B
                        if myAllele != "B":
                            col_data[rowSNP] = "B"
                            NbSNP_Imputed += 1
                    elif freqA < maxFreqAinH and freqB < maxFreqBinH:
                        if freqH >= minFreqH:#'it's a HTZ block
                            if myAllele != "H":
                                col_data[rowSNP] = "H"
                                NbSNP_Imputed += 1
    return(col_data,NbSNP_Imputed)
    
    
def imputation_met3(col_data,WindowSize):
    NbSNP_Imputed=0
    
    for myRow in range(WindowSize):
          #'compare SNP in the middle of the window with the window major allele  
        rowSNP = myRow
        myAllele = col_data[rowSNP]
        myRow1 = rowSNP - WindowSize #'start of window
        myRow2 = rowSNP + WindowSize  #'end of window
        #'handle chr ends
        if myRow1 < 1:
            myRow1 = 1
            myRow2 = rowSNP + (rowSNP - 1) #'to avoid unbalanced windows (otherwise it suppresses small terminal chunks)
        if myRow2 > len(col_data):
            myRow2 = len(col_data)
            myRow1 = rowSNP - (len(col_data) - rowSNP)
                    #'calculate frequencies in window
        nbA = 0; nbB = 0; nbH = 0; nbMD = 0; Nwindow = 0
        for thisRow in range(myRow1,myRow2):
            if col_data[rowSNP] == "A":
                nbA += 1
            elif col_data[rowSNP] == "B":
                nbB += 1
            elif col_data[rowSNP] == "H":
                nbH += 1
            elif col_data[rowSNP] == "-":
                        nbMD += 1
            
        Nwindow = nbA + nbB + nbH + nbMD
        if Nwindow == 0: #this can happen when many consecutive missing data
            col_data[rowSNP] = myAllele
        else:
            freqA = nbA / Nwindow; freqB = nbB / Nwindow; freqH = nbH / Nwindow
            #'compare SNP in the middle of the window with the window major allele
            col_data[rowSNP] = myAllele
                    
            if freqA > freqB and  freqA > freqH:#'it should be an A
                if myAllele != "A":
                    col_data[rowSNP] = "A"
                    NbSNP_Imputed += 1
                elif freqB > freqA and  freqB > freqH:#'it should be an B
                    if myAllele != "B":
                        col_data[rowSNP] = "B"
                        NbSNP_Imputed += 1
                elif freqH > freqA and freqH > freqB:
                    if myAllele != "H":
                        col_data[rowSNP] = "H"
                        NbSNP_Imputed += 1
                    
                    
    for myRow in range(WindowSize):
    
        #chr start
        myAllele = col_data[myRow]
        
        for myRow2 in range(myRow + 1 ,len(col_data)):
            myAllele2 = col_data[myRow2]
            if myAllele2 != "-":
                col_data[myRow]=myAllele2
                NbSNP_Imputed = NbSNP_Imputed + 1
                
    for myRow in range(len(col_data)-1,WindowSize,-1):
        myAllele = col_data[myRow]
        if myAllele == "-":
            for myRow2 in range(myRow -1 ,0,-1):
                myAllele2 = col_data[myRow2]
            if myAllele2 != "-":
                col_data[myRow]=myAllele2
                NbSNP_Imputed = NbSNP_Imputed + 1    
                
                
    return(col_data,NbSNP_Imputed)                
