# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 00:11:55 2020

@author: cami_

"""
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm


def imp_stats(step,n_chr,chr_dir,char_sep,out_dir):
    
    plt.ioff()
    
    step2=step
    steps2 = list(range(2,10+1)) 
    for step2 in steps2:
        if not step in steps2:
            raise ValueError("Please use step 2 to 10",)
     
        
    if 1 < step <=7:
        chr_dir3 = chr_dir+"/workspace"+"/"+"prefiltering"
    else:
        chr_dir3 = chr_dir+"/workspace"+"/"+"post_imputation"
  
    
    
    export_file_path_filt = (chr_dir3+"/"+"Chr"+str(n_chr)+"_step_"+str(step)+".txt")
    graph_file = (out_dir+"/"+"PLOT_LINE_"+"Chr"+str(n_chr)+"_step_"+str(step)+".pdf")
    _summary_list = []
    split_file= open(export_file_path_filt)
    for line in split_file:#tqdm(split_file):
         
          
          
          
        _data =  split_file.readline().strip().split()
        if len(_data)>0:
            _data_f = [{'Position': _data[0].split(char_sep)[1],
                        'MD': _data.count("-")/len(_data)*100,
                        'HTZ': _data.count("H")/len(_data)*100, 
                        'A':_data.count("A")/len(_data)*100, 
                        'B':_data.count("B")/len(_data)*100}] 

            _data_f = pd.DataFrame(_data_f)   
            _summary_list.append(_data_f)
        
    _summary_file = pd.concat(_summary_list)
    _summary_file.index = _summary_file["Position"]
    _summary_file=_summary_file.drop(columns='Position',)
    
    
   
    plt.figure(figsize=(150, 30))
    #plt.figure()
    plt.plot(lines =_summary_file.plot.line())
    plt.title('Chromosome'+" "+str(n_chr)+" "+"step"+" "+str(step), fontsize=20)
    plt.xticks(rotation=90)
    plt.ylabel('(%)', fontsize=18)
    plt.xlabel('Chromosome position', fontsize=18)
    plt.legend(bbox_to_anchor=(.8,-.4),ncol=5)
    #plt.show()
    plt.savefig(graph_file,dpi=300,optimize =True, papertype="legal",bbox_inches="tight") 
    plt.close()   
    return("plotted")
 

n_chrs = list(range(1,12+1))
n_chr=1
steps=list(range(2,10))
char_sep="_"
out_dir = "E:"
chr_dir=chr_dir

for i in tqdm(range(len(n_chrs)),desc="Plotting chromosome"):
    for j in tqdm(range(len(steps)),desc="Plotting all steps"):
        n_chr=n_chrs[i]
        step=steps[j]
        imp_stats(step,n_chr,chr_dir,char_sep,out_dir) #step 9

#x= imp_stats(step,n_chr,chr_dir,char_sep,out_dir)