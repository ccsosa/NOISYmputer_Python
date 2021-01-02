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
import copy
from tqdm import tqdm
import more_itertools as mit
import csv


def small_chunks_step_chr(chr_dir2,n_chrs,smallChunkMaxSize,minEnvirRecRate,chunkEnvironmentMultiplier,maxChunkRecRate,nbOfParents,popSize,mappingFunction,chr_sep,popType,log_dir):
    for i in tqdm(range(len(n_chrs)),desc="Small chunks"):
        n_chr=str(n_chrs[i])
        small_chunks_step(chr_dir2,n_chr,smallChunkMaxSize,minEnvirRecRate,chunkEnvironmentMultiplier,maxChunkRecRate,nbOfParents,popSize,mappingFunction,chr_sep,popType,log_dir)
        
    return("DONE!")

def small_chunks_step(chr_dir2,n_chr,smallChunkMaxSize,minEnvirRecRate,chunkEnvironmentMultiplier,maxChunkRecRate,nbOfParents,popSize,mappingFunction,chr_sep,popType,log_dir):    
   mapping_List = ["Kosambi","Haldane","none"] 
   for map_f in mapping_List:
        if not mappingFunction in mapping_List:
            raise ValueError("No valid mapping function used. Choose a valid option (Kosambi,Haldane,none)")
    
   export_file_path_filt = (chr_dir2+"/"+"Chr"+n_chr+"_step_9.txt")
   export_file_path_filt2 = (chr_dir2+"/"+"Chr"+n_chr+"_step_10.txt")
       
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
   for i in prange(split_file.shape[1]):
       #print(i)
       x = Chunk_col(split_file,chr_pos_num,smallChunkMaxSize,chunkEnvironmentMultiplier,minEnvirRecRate,maxChunkRecRate,nbOfParents,mappingFunction,popSize,popType,i)
       split_file2.iloc[:,(i+1)] = x
   split_file2["#"] = chr_pos
   split_file2.to_csv(export_file_path_filt2,index = False, header=True,sep=" ")
   split_file = pd.read_csv(export_file_path_filt, sep=" ")
   scl = changes_count(split_file2,split_file)
   data =[["step","Small chunks"],
           ["chr",n_chr],
           ["changes",scl]]
   log_file_df = pd.DataFrame(data,columns=["item","status"])
   log_file_df.to_csv(log_dir+"/"+"Chr"+n_chr+"_step9.log",index = False, header=True)

def freq_chunk_preprocess(split_file):
    freqCounts = np.zeros((split_file.shape[0], 4))
    for i in range(split_file.shape[0]):
        freqCounts[i,0] = np.sum(split_file[i,:]=="A")
        freqCounts[i,1] = np.sum(split_file[i,:]=="B")
        freqCounts[i,2] = np.sum(split_file[i,:]=="H")
        freqCounts[i,3] = np.sum(split_file[i,:]=="-")

    return(freqCounts)

freqCounts  = freq_chunk_preprocess(split_file=split_file)


    
def col_chunk_preprocess(split_file):
    colCounts = []
    for i in range(split_file.shape[1]):
        x = list(split_file[:,i])
        x2 = np.where(np.asarray(x)=="-")[0]
        x2 = [item.tolist() for item in x2]
        x2 = [list(group) for group in mit.consecutive_groups(x2)]
        
        for j in range(len(x2)):
            colCounts_df = np.array(
            [i,
            np.min(x2[j]),
            np.max(x2[j]),
            len(x2[j])
            ])
            colCounts.append(colCounts_df)
    #colCounts = np.concatenate(colCounts)
    colCounts = np.stack(colCounts, axis=0)

    return(colCounts)

colCounts  = col_chunk_preprocess(split_file=split_file)
#ind = np.lexsort((colCounts[:,0], colCounts[:,3]))
#colCounts = colCounts[ind]
np.savetxt("D:/counts.csv", colCounts, delimiter=",")

            
def Chunk_col(split_file,chr_pos_num,smallChunkMaxSize,chunkEnvironmentMultiplier,minEnvirRecRate,maxChunkRecRate,nbOfParents,mappingFunction,popSize,popType,i):
   # i = 0
    col_data = split_file[:,i]
    col_data = copy.deepcopy(col_data)
    
    _ignore_cols = copy.deepcopy(nbOfParents)
    _ignore_cols.append(i)
    _ignore_cols=list(np.unique(_ignore_cols))
    
    _col_sub0=[i for i in range(split_file.shape[1]) if i not in _ignore_cols]
    _col_sub0 = np.array(copy.deepcopy(split_file[:,_col_sub0]))
                             
                                   
    #myCol=i
    for myRow in range(1,len(col_data)-2):
        print("row= "+str(myRow))

        myAllele = col_data[myRow]
        myAllele1 = col_data[myRow-1] #allele just before the candidate chunk
        
        if myAllele==myAllele1:
            continue
        else:
            #myRow+=1
            #pass
        #else:#possible start of a chunk detected
            #print("row= "+str(myRow)+ "| " )
            myStartRowChunk=myRow #store start position
            myStartPosChunk=chr_pos_num[myRow] #in BP
            #myRow2 = myRow+1
            for myRow2 in range(myRow+1,len(col_data)-1):
            
                #print("row= "+str(myRow)+ "| " +"myRow2 : "+str(myRow2))
                #print(myRow2)
                
                if myRow2> len(col_data)-1:
                    break
                else:
                    myAllele2 =  col_data[myRow2] #allele of SNP just after end of candidate chunk
                    if myAllele2!=myAllele1: #possible end of chunk detected
                        myStopRowChunk =myRow2-1
                        thisChunkSize = myStopRowChunk - myStartRowChunk + 1 # chunk size in number of rows
                        chunkEnvironmentSize = int(thisChunkSize * chunkEnvironmentMultiplier) + 1
                        if thisChunkSize <= smallChunkMaxSize: #it's a chunk of acceptable size
                            #myChunkCol = i
                            myStopPosChunk = chr_pos_num[myStopRowChunk] #in BP
                            chunkSizeBP = myStopPosChunk - myStartPosChunk
                            EnvironmentHomogeneous = True
                            startEnvirRow = myStartRowChunk - chunkEnvironmentSize
                            stopEnvirRow = myStopRowChunk + chunkEnvironmentSize
                            
                            if myAllele2 == myAllele1 or myAllele2 == "-" or myAllele1 == "-": #verify that the environment is homogeneous
                                if myAllele1 !="-": 
                                    envirAllele = myAllele1
                                if myAllele2 !="-": 
                                    envirAllele = myAllele2
                                #myRow3 = myStartRowChunk-1 
                                for myRow3 in range(myStartRowChunk-1,startEnvirRow,-1):#verify homogeneity up
                                    #print("row= "+str(myRow)+ "| " +"myRow3 1: "+str(myRow3))
                                    if myRow3 < 1: 
                                        EnvironmentHomogeneous = False
                                        break
                                    if col_data[myRow3]!=envirAllele:
                                        if  col_data[myRow3]!="-":   #some MD can be around
                                            EnvironmentHomogeneous = False
                                    
                                    
                                #myRow3 = myStopRowChunk+1   
                                for myRow3 in range(myStopRowChunk+1,stopEnvirRow): #verify homogeneity down
                                    #print("row= "+str(myRow)+ "| " +"myRow3 2: "+str(myRow3))
                                    if myRow3 >len(col_data): 
                                        EnvironmentHomogeneous = False
                                        break
                                    if  col_data[myRow3]!=envirAllele: #not homogeneous, stop
                                        if  col_data[myRow3]!="-":   
                                            EnvironmentHomogeneous = False
                                    #myRow3+= 1
                               
                                if EnvironmentHomogeneous==False:
                                    continue
                                else:# EnvironmentHomogeneous==True:
                                   
                                    """
                                    #'now we need the local probability of recombination to decide if we replace alleles
                                    start_E = math.ceil(startEnvirRow)
                                    stop_E = math.ceil(stopEnvirRow)
                                    if math.ceil(start_E)<0:
                                        startEnvirRow = 1
                                    else:
                                        startEnvirRow = start_E
                                        if math.ceil(stop_E)<0:
                                                stopEnvirRow = 1
                                        else:
                                            stopEnvirRow= stop_E
                                    """
                                 
                                    if np.sign(startEnvirRow)==-1:
                                        startEnvirRow = 0 
                                    else:
                                        startEnvirRow = startEnvirRow 
                                    
                                    startEnvirPos = chr_pos_num[startEnvirRow]
                                    stopEnvirPos = chr_pos_num[stopEnvirRow]
                                    EnvirSizeBP = stopEnvirPos - startEnvirPos #envir size in bp
                                
                                
                                    #method 2, calculate local recombination rate in the chunk environment
                                    EnvircMSize = 0
                                    for myRow3 in range(startEnvirRow,stopEnvirRow-1):
                                        #print("row= "+str(myRow)+ "| " +"myRow3 3: "+str(myRow3))
                                        myCol2= i
                                        recFrac = 0
                                        F1RecFreq = 0
                                        F1cM = 0
                                        _nbAA = 0
                                        _nbAB = 0
                                        _nbAH = 0
                                        _nbBA = 0
                                        _nbBB = 0
                                        _nbBH = 0
                                        _nbMD = 0
                                        _nbHA = 0
                                        _nbHB = 0
                                        _nbHH = 0
                                        _Npop = 0
                                          
                                       
                                        for myCol2 in range(_col_sub0.shape[1]-1):
                                            #print("row= "+str(myRow)+ "| " +"myCol2 1: "+str(myCol2))
                                            __myAlleleSNP1 = _col_sub0[myRow3,myCol2]
                                            __myAlleleSNP2 = _col_sub0[int(myRow3 +1),myCol2]
                                            if __myAlleleSNP1 ==  "H":
                                                if __myAlleleSNP2 == "H":
                                                    _nbHH += 1
                                                elif __myAlleleSNP2 == "A":
                                                    _nbHA +=1
                                                elif __myAlleleSNP2 == "B":
                                                    _nbHB +=1
                                               
                                            elif __myAlleleSNP1 == "A":
                                                if __myAlleleSNP2=="A":
                                                   _nbAA +=1
                                                elif __myAlleleSNP2=="H":
                                                    _nbAH+=1
                                                elif __myAlleleSNP2=="B":
                                                    _nbAB +=1
                                                
                                            elif __myAlleleSNP1=="B":
                                                if __myAlleleSNP2=="B":
                                                    _nbBB +=1
                                                elif __myAlleleSNP2=="H":
                                                    _nbBH=+1
                                                elif __myAlleleSNP2=="A":
                                                    _nbBA +=1
                                                elif __myAlleleSNP1=="-":
                                                    _nbMD+=1
                                            
                                        if popType == "SSD": #'population size for the current SNP pair
                                            _Npop = _nbAA + _nbBB + _nbAB + _nbBA + _nbMD
                                        elif popType == "BC1":
                                            _Npop = _nbAA + _nbAH + _nbHH + _nbHA + _nbMD
                                        elif popType == "DH":
                                            _Npop = _nbAA + _nbAB + _nbBB + _nbBA + _nbMD
                                        elif popType == "F2":
                                            _Npop = _nbAA + _nbBB + _nbHH + _nbAB + _nbBA + _nbAH + _nbBH + _nbHA + _nbHB + _nbMD


                                        if _Npop > 0: # 'recombination fraction
                                            if popType == "SSD":
                                                recFrac = (_nbAB + _nbBA) / _Npop  #'for SSD, we could use the htz info like in Calculate_final_map 'this is the rec frac observed in SSD, not converted to F1 value
                                                RecFreq= CountBKPT(split_file)/popSize
                                                try:
                                                    F1RecFreq = float(RecFreq) / float(2 - 2 * recFrac)
                                                except ZeroDivisionError:
                                                    F1RecFreq = 0
                                                    #F1RecFreq = RecFreq / (2 - 2 * recFrac) #'we could use Martin & Hospital correction but that wouldn't change much
                                            elif popType == "F2": # 'EM algorithm
                                                F1RecFreq = RecFreq_F2(_nbAA, _nbBB, _nbHH, _nbAB, _nbBA, _nbAH, _nbBH, _nbHA, _nbHB, _Npop)
                                            elif popType== "BC1": #Then
                                                F1RecFreq = (_nbAH + _nbHA) / _Npop
                                            elif popType == "DH" :
                                                F1RecFreq = (_nbAB + _nbBA) / _Npop
                                        else:
                                            F1RecFreq = 0
                                
                                        if mappingFunction == "Kosambi": #'conversion to centimorgans
                                            F1cM = Kosambi(F1RecFreq)
                                        elif mappingFunction == "Haldane":
                                            F1cM = Haldane(F1RecFreq)
                                        elif mappingFunction == "none" :
                                            F1cM = 100 * F1RecFreq    
                                        EnvircMSize = EnvircMSize + F1cM
                                        
                                    envirRecRate = EnvircMSize / (EnvirSizeBP / 1000000) #'cM/mega bp
                            
                                    #'we also check that the chunk itself isn't too high in recombination, independent to its environment
                                    #                'we have chunkSizeBP     
                            
                                    ThisChunkcMSize = 0
                                    for myRow3 in range(myStartRowChunk - 1, myStopRowChunk):
                                        #print("row= "+str(myRow)+ "| " +"myRow3 4: "+str(myRow3))
                                        recFrac = 0
                                        F1RecFreq = 0
                                        F1cM = 0
                                        _nbAA = 0
                                        _nbAB = 0
                                        _nbAH = 0
                                        _nbBA = 0
                                        _nbBB = 0
                                        _nbBH = 0
                                        _nbMD = 0
                                        _nbHA = 0
                                        _nbHB = 0
                                        _nbHH = 0
                                        _Npop = 0
                                        
   
                                        
                                        #col_sub_2 =np.array([split_file[j,m] for m in range(split_file.shape[1]) if m not in nbOfParents])
                                        for myCol2 in range(_col_sub0.shape[1]-1):
                                            #print("row= "+str(myRow)+ "| " +"myCol2 2: "+str(myCol2))
                                            
                                            _myAlleleSNP1 = _col_sub0[myRow3,myCol2]
                                            _myAlleleSNP2 = _col_sub0[int(myRow3 + 1),myCol2]
                                            if _myAlleleSNP1 ==  "H":
                                                if _myAlleleSNP2 == "H":
                                                    _nbHH += 1
                                                elif _myAlleleSNP2 == "A":
                                                    _nbHA +=1
                                                elif _myAlleleSNP2 == "B":
                                                    _nbHB +=1
                                               
                                            elif _myAlleleSNP1 == "A":
                                                if _myAlleleSNP2=="A":
                                                    _nbAA +=1
                                                elif _myAlleleSNP2=="H":
                                                    _nbAH+=1
                                                elif _myAlleleSNP2=="B":
                                                    _nbAB +=1
                                            elif _myAlleleSNP1=="B":
                                                if _myAlleleSNP2=="B":
                                                    _nbBB +=1
                                                elif _myAlleleSNP2=="H":
                                                    _nbBH=+1
                                                elif _myAlleleSNP2=="A":
                                                    _nbBA +=1
                                                elif _myAlleleSNP1=="-":
                                                    _nbMD+=1
                                            
                                        if popType == "SSD": #'population size for the current SNP pair
                                            _Npop = _nbAA + _nbBB + _nbAB + _nbBA + _nbMD
                                        elif popType == "BC1":
                                            _Npop = _nbAA + _nbAH + _nbHH + _nbHA + _nbMD
                                        elif popType == "DH":
                                            _Npop = _nbAA + _nbAB + _nbBB + _nbBA + _nbMD
                                        elif popType == "F2":
                                           _Npop = _nbAA +_nbBB + _nbHH + _nbAB + _nbBA + _nbAH + _nbBH + _nbHA + _nbHB + _nbMD


                                        if _Npop > 0: # 'recombination fraction
                                            if popType == "SSD":
                                                recFrac = float((_nbAB + _nbBA) / _Npop) # 'for SSD, we could use the htz info like in Calculate_final_map 'this is the rec frac observed in SSD, not converted to F1 value
                                                RecFreq = CountBKPT(split_file)/split_file.shape[1]
                                                try:
                                                    F1RecFreq = float(RecFreq) / float(2 - 2 * recFrac)
                                                except ZeroDivisionError:
                                                    F1RecFreq = 0
                                            #'we could use Martin & Hospital correction but that wouldn't change much
                                            elif popType == "F2": # 'EM algorithm
                                                F1RecFreq = RecFreq_F2(_nbAA, _nbBB, _nbHH, _nbAB, _nbBA, _nbAH, _nbBH, _nbHA, _nbHB, _Npop)
                                            elif popType == "BC1": #Then
                                                F1RecFreq = (_nbAH + _nbHA) / _Npop
                                            elif popType == "DH" :
                                                F1RecFreq = (_nbAB + _nbBA) / _Npop
                                            else:
                                                F1RecFreq = 0
                                
                                        if mappingFunction == "Kosambi": #'conversion to centimorgans
                                            F1cM = Kosambi(F1RecFreq)
                                        elif mappingFunction == "Haldane":
                                            F1cM = Haldane(F1RecFreq)
                                        elif mappingFunction == "none" :
                                            F1cM = 100 * F1RecFreq    
                                        
                                        ThisChunkcMSize = ThisChunkcMSize + F1cM
                                     
                                    if chunkSizeBP > 0:
                                        ThisChunkRecRate = ThisChunkcMSize / (chunkSizeBP / 1000000) #'cM/mega bp
                                    else:
                                         ThisChunkRecRate = maxChunkRecRate + 1
                            
                                    if envirRecRate <= minEnvirRecRate or ThisChunkRecRate >= maxChunkRecRate: #'not enough local recomb, or chunk rec rate too high
                                        #Debug.Print myCol, myStartPosChunk, myStopPosChunk, envirRecRate, EnvirSizeBP, ThisChunkRecRate, chunkSizeBP
                                        #           'proceed to allele replacement
                                        for myRow3 in range(myStartRowChunk,myStopRowChunk): #'replace the chunk allele by its environment allele
                                            col_data[myRow3] = myAllele1
                                
    return col_data 



def CountBKPT(split_file):
    for x in range(split_file.shape[0]):
        col_sub_BKPT = split_file[x,:]
        CountBKPT=0
        for y in range(len(col_sub_BKPT)-1):
            myRow =y 
            myRow2 = y + 1
            myAllele = col_sub_BKPT[myRow]
            myAllele2 = col_sub_BKPT[myRow2]
            if myAllele != myAllele2 and myAllele != "-" and myAllele2 != "-": #'bkpt
                CountBKPT = CountBKPT + 1
    return(CountBKPT)

def RecFreq_F2(nbAA, nbBB, nbHH, nbAB, nbBA, nbAH, nbBH, nbHA, nbHB, Npop):
    if Npop > 0:
        invNpop = 1 / Npop
        p = (2 * (nbAB + nbBA) + (nbAH + nbHA + nbHB + nbBH) / 2) * invNpop
        #'NbIter = 0
        Adjust = 1
        while Adjust > 0.0001:
            t = (p ** 2) / ((1 - p) ** 2 + p ** 2)
            pNew = (2 * (nbAB + nbBA) + (nbAH + nbHA + nbHB + nbBH) + 2 * t * nbHH) * invNpop * 0.5
            #'NbIter = NbIter + 1
            Adjust = abs((pNew - p))
            p = pNew
        next
        RecFreq_F2 = pNew
    else:
        RecFreq_F2 = 0
    return(RecFreq_F2)


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


