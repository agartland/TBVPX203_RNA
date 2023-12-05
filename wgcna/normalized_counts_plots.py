import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
import itertools
import os
from os.path import join as opj
import sys
import os
from pngpdf import PngPdfPages
from biplot import biplot
from myboxplot import swarmbox
import matplotlib.pyplot as plt

# Your plotting code here

# Other code or plot customization

sns.set_style('whitegrid')


normalized_data_all = pd.read_csv('X:/fast/gilbert_p/hsrivast/TBVPX_203/updated_wgcna_pca_eigen/normalized_data_all.csv', index_col=0)


WGCNA_modules_gene = pd.read_csv('X:/fast/gilbert_p/hsrivast/TBVPX_203/updated_wgcna_pca_eigen/WGCNA_modules_stricter_DEG_expression.csv')
normalized_modules = pd.DataFrame()

rx = pd.read_csv('X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/trt_pubid_2022-DEC-19.csv')
mod_wgcna=["turquoise","blue","green","yellow","black","red","brown"]


for i in mod_wgcna:
    first_mod = WGCNA_modules_gene[i].dropna()
    subset_df = normalized_data_all.loc[first_mod]
    normalized_modules[i]=subset_df.mean()
 
    
eg=normalized_modules
eg['sampleid']=eg.index
eg = eg.assign(Visit=eg['sampleid'].map(lambda s: f"D{s.split('_')[-1]}"),
               ptid=eg['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])))

eg = pd.merge(eg, rx[['Treatment_Group','ptid', 'Treatment']], how='left', on='ptid')



modules = ['turquoise','blue','green','yellow','black','red','brown']




with PngPdfPages('normalized_count_2.pdf') as pdf:
    sns.clustermap(eg[modules].corr(), metric='correlation', annot=True)
    pdf.savefig(plt.gcf())
    plt.close(plt.gcf())

    sns.clustermap(eg[modules], metric='correlation')
    pdf.savefig(plt.gcf())
    plt.close(plt.gcf())

    for trt in [('2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)', '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)')]:
        for m in modules:
            plotdf = eg.loc[eg['Treatment_Group'].isin(trt)]
            figh = plt.figure(figsize=(10,10))
            swarmbox(x='Visit', y=m, connect=True, connect_on=['ptid'], data=plotdf, order=['D0', 'D3', 'D56', 'D59', 'D63', 'D112', 'D168'])
            plt.title(f'{trt[0]}\n{trt[1]}')
            plt.xlabel('Visits')
            plt.ylabel(f'{m}_module_score')  
            pdf.savefig(figh)

    
            