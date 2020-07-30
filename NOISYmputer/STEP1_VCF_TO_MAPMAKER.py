# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:02:16 2020
NOISYmputer python step 0 VCF to mapmaker formats
@author: Chrystian Sosa
@Correspondence author: Dr Mathias Lorieux
"""

import pandas as pd
import platform
import subprocess
import os.path

def transform_vcf_to_mapmaker(path_to_addon,popType,vcf_file,parent1_id,parent2_id,output_ABH_file):

#IMPORTING MODULES
    

#GATHERING OPERATIVE SYSTEM
    op_sys = platform.system()
#CHECKING IF OUTPUT FILE WAS MADE PREVIOUSLY
    if os.path.isfile(output_ABH_file):
        print ("VCF file converted to mapmaker format using "+popType+" population "+ " was made previously, loading...")
        output = pd.read_csv(output_ABH_file, sep=" ", skiprows=2)
        return(output)
#RUNNING JAVA SHELL COMMMAND USING SUBPROCESS
    else: 
        print("Running MapDistoAddonsMT_v5.jar")
        process = subprocess.run('java -jar '
                          + path_to_addon +
                          ' -T ImportVCF -OS '
                          + op_sys.lower()+ 
                          ' -d '
                          + popType + 
                          ' -f ' + 
                          vcf_file +
                          ' -p '+
                          parent1_id +','+parent2_id +
                          ' -o '+output_ABH_file,
                          check=True,
                          shell=True,
                          stdout=subprocess.PIPE,
                          universal_newlines=True
                          )
#LOADING FINAL RESULTS CHECKING IF THE TRANSFORMED MATRIX WAS DONE CORRECTLY AND LOADING INTO A PANDAS DATAFRAME
        if os.path.isfile(output_ABH_file):
            print ("VCF file converted to mapmaker format using "+popType+" population")
            output = pd.read_csv(output_ABH_file, sep=" ", skiprows=2)
            return(output)
        else:
#GETTING ERROR MESSAGE WHEN THERE IS NOT TRANSFORMED MATRIX IN output_ABH_file PATH
            print ("File not exist, try again")
            output = print("error!, VCF to mapmaker format failed, check parameters and try again")
            return(output)
        
#RUNNING TEST
            
#path_to_addon = "C:/Users/cami_/Documents/MapDisto_workspace/MapDisto_plugins/MapDistoAddonsMT_v5.jar"
#popType = "ri_self"
#vcf_file ="C:/Users/cami_/Documents/MapDisto_data/Rice_GBS.vcf"
#parent1_id = "ID152bH10-P2_CGTACG-GTGGCC"
#parent2_id = "ID152bH11-P2_GAGTGG-GTGGCC"
#output_ABH_file = "E:/rice_gbs2.txt"
#out = transform_vcf_to_mapmaker(path_to_addon,popType,vcf_file,parent1_id,parent2_id,output_ABH_file)
#chr_dir = "E:/CHR"

