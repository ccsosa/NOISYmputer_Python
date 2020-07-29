# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 21:30:18 2020

@author: cami_
"""

import os
import sys
sys.path.append("E:/CHR/TEST/SCRIPTS")

#PREFILTERING
from  STEP1_VCF_TO_MAPMAKER import transform_vcf_to_mapmaker
from  STEP2_SPLIT_CHROMOSOME import split_by_chromosome
from  STEP3_FILTERING_SNP_STEP import filtering_snp_step_chr
from  STEP3_FILTERING_SNP_STEP import filtering_snp_step
from  STEP4_FILTERING_MISSING_DUPS import remove_dup_snp_step_chr
from  STEP4_FILTERING_MISSING_DUPS import remove_dup_snp_step
from  STEP5_CORRECT_ERRORS import correct_chr_read_chr
from  STEP5_CORRECT_ERRORS import correct_chr_read
from  STEP5_CORRECT_ERRORS import correct_loci_col
from  STEP6_FILLING_MISSING_DATA import chr_read_missing_data_chr
from  STEP6_FILLING_MISSING_DATA import chr_read_missing_data
from  STEP6_FILLING_MISSING_DATA import correct_loci_col_miss
from  STEP7_MARK_INCOHERENT import correct_loci_row_chr
from  STEP7_MARK_INCOHERENT import correct_loci_row
#POST-IMPUTATION
from  STEP8_IMPUTATION import imputation_chr
from  STEP8_IMPUTATION import imputation
from  STEP8_IMPUTATION import imputation_col
from  STEP8_IMPUTATION import imputation_0
from  STEP8_IMPUTATION import imputation_met1
from  STEP8_IMPUTATION import imputation_met2
from  STEP8_IMPUTATION import imputation_met3
from  STEP9_SINGLETONS_BREAKS import correct_singletons_break_chr
from  STEP9_SINGLETONS_BREAKS import correct_singletons_break
from  STEP9_SINGLETONS_BREAKS import mark_single_rows
from  STEP10_SMALL_CHUNKS import small_chunks_step_chr
from  STEP10_SMALL_CHUNKS import small_chunks_step
from  STEP10_SMALL_CHUNKS import Chunk_col
from  STEP10_SMALL_CHUNKS import CountBKPT
from  STEP10_SMALL_CHUNKS import RecFreq_F2
from  STEP10_SMALL_CHUNKS import Kosambi
from  STEP10_SMALL_CHUNKS import Haldane
#AUXILIAR EXPRESSIONS
from  STEP0_COUNT_CHANGES import changes_count
from  STEP0_PLOT_CHR import chr_plot_all
from  STEP0_PLOT_CHR import chr_plot
from  STEP0_IMP_STATS import imp_stats



def NOISYmputer_pipeline(pop_markers_params,prefiltering_SNPs_params,filtering_incoherent_SNPs_params,imputation_params_step1,imputation_params_step2,imputation_params_step3,Breakpoint_singletons_params,improbable_chunks_params_params,Alien_segments_params):
    
    
    n_chrs=  pop_markers_params[2][1]
    path_to_addon=pop_markers_params[3][1]
    vcf_file=pop_markers_params[4][1]
    chr_dir=pop_markers_params[5][1]
    parent1_id=pop_markers_params[6][1]
    parent2_id=pop_markers_params[7][1]
    output_ABH_file=pop_markers_params[7][1]
    
   
        
    #DEFINING AND CHECKING POP TYPE
    popType= pop_markers_params[0][1]

    popTypeList = ["SSD","BC1","DH","F2"] 
    for map_p in popTypeList:
        if not popType in popTypeList:
            raise ValueError("No valid population type. Choose a valid option (SSD,BC1,DH,F2)")
  
    
    #DEFINING AND CHECKING MAPPING FUNCTION
    mappingFunction= pop_markers_params[11][1]
    mapping_List = ["Kosambi","Haldane","none"] 
    for map_f in mapping_List:
        if not mappingFunction in mapping_List:
            raise ValueError("No valid mapping function used. Choose a valid option (Kosambi,Haldane,none)")
  
    
    workspace_dir = chr_dir+"/workspace"
    if(os.path.isdir(workspace_dir)):
        print(workspace_dir, "created")
    else:
            os.mkdir(workspace_dir)
            
    log_dir = workspace_dir+"/log"
    if(os.path.isdir(log_dir)):
        print(log_dir, "created")
    else:
            os.mkdir(log_dir)

    chr_dir1=chr_dir+"/workspace"+"/"+"prefiltering"
    #CHECKING IF CHROMOSOME FOLDER EXISTS
    if(os.path.isdir(chr_dir1)):
        print(chr_dir1, "created")
    else:
            os.mkdir(chr_dir1)
    
    
    chr_dir2=chr_dir+"/workspace"+"/"+"post_imputation"
    #CHECKING IF CHROMOSOME FOLDER EXISTS
    if(os.path.isdir(chr_dir2)):
        print(chr_dir1, "created")
    else:
            os.mkdir(chr_dir2)
            
            
    ######################################################################################        
    #PREFILTERING
    #STEP 1
    x = transform_vcf_to_mapmaker(path_to_addon,popType,vcf_file,parent1_id,parent2_id,output_ABH_file)
    x= split_by_chromosome(chr_dir1,input=x)
    
    #STEP 3 Pre-filtering SNPs 1

     
    maxFreqMD=prefiltering_SNPs_params[0][1]
    maxFreqH=prefiltering_SNPs_params[1][1]
    minFreqA=prefiltering_SNPs_params[2][1]
    minFreqB=prefiltering_SNPs_params[3][1]
    
    x = filtering_snp_step_chr(chr_dir1,n_chrs,maxFreqMD,maxFreqH,minFreqA,minFreqB,log_dir)
    #STEP 4 Pre-filtering SNPs 2
    x = remove_dup_snp_step_chr(chr_dir1,n_chrs,"missing",log_dir)    
    #STEP 5 Filtering incoherent SNPs
    WindowSize=7
    x = correct_chr_read_chr(chr_dir1,n_chrs,popType,WindowSize,log_dir)
    #STEP 6 filling missing data
    x = chr_read_missing_data_chr(chr_dir1,[n_chrs[0]],log_dir)
    #STEP 7
    WindowSize=filtering_incoherent_SNPs_params[0][1]    
    Chi2threshold = filtering_incoherent_SNPs_params[1][1]    
    x= correct_loci_row_chr(chr_dir1,n_chrs,Chi2threshold,WindowSize,log_dir)
    
    ######################################################################################
    #IMPUTATION
    #STEP 8
    chr_sep=pop_markers_params[9][1]
    x=imputation_chr(chr_dir1,chr_sep,n_chrs,chr_dir2,imputation_params_options,imputation_params_step1_params,imputation_params_step2_params,imputation_params_step3_params,log_dir)
    ######################################################################################
    #POSTIMPUTATION
    #STEP 9
    chr_sep= pop_markers_params[9][1]
    x=correct_singletons_break_chr(n_chrs,chr_dir2,chr_sep,log_dir)
    
    #STEP 10
    
    

    popSize=pop_markers_params[1][1]
    mappingFunction=pop_markers_params[11][1]
    chr_sep= pop_markers_params[9][1]
    popType=pop_markers_params[0][1]
    smallChunkMaxSize =improbable_chunks_params[0][1]
    chunkEnvironmentMultiplier = improbable_chunks_params[1][1]
    minEnvirRecRate =improbable_chunks_params[1][1]
    maxChunkRecRate = improbable_chunks_params[3][1]
    nbOfParents= pop_markers_params[12][1]
    x=small_chunks_step_chr(chr_dir2,n_chrs,smallChunkMaxSize,minEnvirRecRate,chunkEnvironmentMultiplier,maxChunkRecRate,nbOfParents,popSize,mappingFunction,chr_sep,popType,log_dir)

    return("DONE!")
    
    
pop_markers_params = [
          ["population type","SSD"],
          ["population size",187],
          ["nb of chromosomes",[1]],#,list(range(1,12+1))],
          ["path_to_addon","C:/Users/cami_/Documents/MapDisto_workspace/MapDisto_plugins/MapDistoAddonsMT_v5.jar"],
          ["vcf_file","C:/Users/cami_/Documents/MapDisto_data/Rice_GBS.vcf"],
          ["chr_dir","E:/CHR/TEST"],
          ["Parent 1 ID in VCF","ID152bH10-P2_CGTACG-GTGGCC"],
          ["Parent 2 ID in VCF","ID152bH11-P2_GAGTGG-GTGGCC"],
          ["output_ABH_file","E:/CHR/rice_gbs2.txt"],
          ["chr-position separator in SNP names","_"],
          ["cM/Mbp",250.0],
          ["mappingFunction","Kosambi"],
          ["nbOfParents",[0,1]]
          ]


prefiltering_SNPs_params= [
          ["max missing data",0.666],
          ["max heterozygosity",0.800],
          ["min f(A)",0.010],
          ["min f(B)",0.010]
          ]

filtering_incoherent_SNPs_params = [
          ["half window size",1000],
          ["local Chi2 threshold",3.84]]



imputation_params_options  = [
          ["method 1",False],
          ["method 2",True],
          ["method 3",False]
          ]


imputation_params_step1_params = [
          ["half window size",50],
          ["shift",1],
          ["min f(A) in A",0.650],
          ["min f(B) in B",0.650],
          ["min f(A) in H",	0.010],
          ["min f(B) in H",	0.010],
          ["Htz transition rate",0.010]
          ]

imputation_params_step2_params = [
          ["half window size",30],
          ["shift",1],
          ["min f(A) in A",0.900],
          ["min f(B) in B",0.900],
          ["min f(A) in H",	0.550],
          ["min f(B) in H",	0.550],
          ["Htz transition rate",0.010]
          ]

imputation_params_step3_params = [
          ["half window size",25],
          ["shift",1]
          ]
        
Breakpoint_singletons_params = [
          ["1/2 window size",15],
          ["shift",1]
          ]

improbable_chunks_params = [
          ["max chunk size (nb of SNPs)",50],
          ["chunk environment multiplier ",0.51],
          ["min local recomb rate in chunk environment (cM/Mbp)",15],
          ["max local recomb rate in chunk (cM/Mbp)",500]
          ]

     =  [
          ["max window size",100],
          ["maxRF",0.01],
          ["slope threshold",10],
          ["min alien segment size",3],
          ["start alien segment offset",2]
          ]


