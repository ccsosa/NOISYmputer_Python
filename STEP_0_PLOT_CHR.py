# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 21:40:59 2020

@author: cami_
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
from tqdm import tqdm


#https://stackoverflow.com/questions/25482876/how-to-add-legend-to-imshow-in-matplotlib

def chr_plot_all(chr_dir,n_chrs,out_dir):
    
    for i in tqdm(range(1,n_chrs+1),desc="..."):
        for j in  [3,4,5,6,7,9]:
            n_chr=str(i)
            step = str(j)
            chr_plot(chr_dir,n_chr,step,out_dir)
    return("DONE!")


def chr_plot(chr_dir,n_chr,step,out_dir):
    
    
    plt.ioff()

    export_file_path_filt2 = (chr_dir+"/"+"Chr"+n_chr+"_step_"+step+".txt")
    graph_file = (out_dir+"/"+"Chr"+n_chr+"_step_"+step+".pdf")
    
    print("plotting"+" | "+"Chromosome:"+str(n_chr)+" | "+ "step:"+ str(step))
    split_file = pd.read_csv(export_file_path_filt2, sep=" ")
    split_file = split_file.drop(columns=['#'])
    split_file = split_file.to_numpy()

    texts=["-","B","H","A"]
    for i in range(split_file.shape[1]):
        #print(str(i))
        split_file[:,i]= np.where(split_file[:,i]=="A",1,split_file[:,i])
        split_file[:,i]= np.where(split_file[:,i]=="H",0,split_file[:,i])
        split_file[:,i]= np.where(split_file[:,i]=="B",-1,split_file[:,i])
        split_file[:,i]= np.where(split_file[:,i]=="-",-2,split_file[:,i])
    split_file=split_file.astype(int)    
    values = np.unique(split_file)

    plt.figure(figsize=(20, 200))
    im=plt.imshow(split_file, interpolation='none', cmap='RdBu_r')
    colors = [ im.cmap(im.norm(value)) for value in values]
    #patches = [ mpatches.Patch(color=colors[i], label="Level {l}".format(l=values[i]) ) for i in range(len(values)) ]
    patches = [ mpatches.Patch(color=colors[i], label="SNP status: {l}".format(l=texts[i]) ) for i in range(len(values)) ]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    #plt.show()
    plt.savefig(graph_file,dpi=300,optimize =True, papertype="legal",bbox_inches="tight") 
    plt.close()   
    return("plotted")


chr_dir = "E:/CHR"  
out_dir = "E:"

""" 
for i in [3,4,5,6,7,9]:
    print(str(i))
    step = str(i)
    chr_plot(chr_dir,n_chr,step,out_dir)
 """  
n_chrs=12
chr_plot_all(chr_dir,n_chrs,out_dir)
    
