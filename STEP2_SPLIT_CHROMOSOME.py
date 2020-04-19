# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:02:16 2020
NOISYmputer python step 0 Splitting mapmaker by choromosomes
@author: Chrystian Sosa
@Correspondence author: Dr Mathias Lorieux
"""
def split_by_chromosome(chr_dir,out):
#CALLING PYTHON MODULES
    import re 
    import os
    import pandas as pd
#CHECKING IF CHROMOSOME FOLDER EXISTS
    if(os.path.isdir(chr_dir)):
        print(chr_dir, "created")
    else:
            os.mkdir(chr_dir)
# SPLITTING DATA BY CHROMOSOME   
    chr_pos = []
    for i in range(len(out)):
        chr_pos.append(re.split('_', out["#"][i])[0])
                                         
    chr_pos = [c.replace("*", "") for c in chr_pos]
    out['chr'] = chr_pos 
    n_chr = out['chr'].unique()
#SPLITTING FILE PER CHROMOSOME
    for i in range(len(n_chr)):
#CHECKING IF CHROMOSOME FILE WAS ALREADY MADE
        export_file_path = (chr_dir+"/"+n_chr[i]+"_step_2.txt")
        if os.path.isfile(export_file_path):
            print(export_file_path+" available")
        else:
#SAVING CHROMOSOME FILE
            x  = out[out.chr == n_chr[i]]
            x.drop(columns=['chr'])
            del x["chr"]
            x.to_csv(export_file_path,index = False, header=True,sep=" ")
    return(print("chromosomes splited!"))

#RUNNING TEST
    
"""
path_to_addon = "C:/Users/cami_/Documents/MapDisto_workspace/MapDisto_plugins/MapDistoAddonsMT_v5.jar"
pop_type = "ri_self"
vcf_file ="C:/Users/cami_/Documents/MapDisto_data/Rice_GBS.vcf"
parent1_id = "ID152bH10-P2_CGTACG-GTGGCC"
parent2_id = "ID152bH11-P2_GAGTGG-GTGGCC"
output_file = "E:/rice_gbs2.txt"
out = transform_vcf_to_mapmaker(path_to_addon,pop_type,vcf_file,parent1_id,parent2_id,output_file)
chr_dir = "E:/CHR"

split_by_chromosome(chr_dir,out)

"""