import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
import itertools
import os
from os.path import join as opj
import sys

from fg_shared import _git, _fg_data

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from biplot import biplot
from myboxplot import swarmbox

sns.set_style('whitegrid')

eg_fn = opj(_fg_data, 'SCRI/TBVPX-203/RNA/2019Dec/Results/eigen_genes.csv')
eg = pd.read_csv(eg_fn).rename({'Unnamed: 0':'sampleid'}, axis=1)

rx_fn = opj(_fg_data, 'SCRI/TBVPX-203/RNA/trt_pubid_2022-DEC-19.csv')
rx = pd.read_csv(rx_fn)

eg = eg.assign(Visit=eg['sampleid'].map(lambda s: f"D{s.split('_')[-1]}"),
               ptid=eg['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])))

eg = pd.merge(eg, rx[['ptid', 'Treatment']], how='left', on='ptid')

modules = ['turquoise', 'brown', 'grey', 'yellow', 'blue']

with PngPdfPages(opj(_fg_data, 'SCRI/TBVPX-203/RNA/2019Dec/Results/eigen_gene_plots.pdf')) as pdf:
    sns.clustermap(eg[modules].corr(), metric='correlation', annot=True)
    pdf.savefig(plt.gcf())
    plt.close(plt.gcf())

    sns.clustermap(eg[modules], metric='correlation')
    pdf.savefig(plt.gcf())
    plt.close(plt.gcf())

    figh = plt.figure(figsize=(10,10))
    for c in modules:
        sns.distplot(eg[c], color=c)
    pdf.savefig(figh)

    figh = plt.figure(figsize=(10, 10))
    axh = figh.add_axes([0.15, 0.15, 0.5, 0.7])
    keep_ind = ~eg['Treatment'].isin([1, 2, 5]) & eg['Visit'].isin(['D56', 'D59', 'D63'])
    biplot(eg.loc[keep_ind, modules], labels=eg.loc[keep_ind, 'Visit'], plotLabels=False, varThresh=0.1)
    pdf.savefig(figh)

    for trt in [(1, 2), (3, 4), (5,)]:
        for m in modules:
            figh = plt.figure(figsize=(15,10))
            plotdf = eg.loc[eg['Treatment'].isin(trt)]
            swarmbox(x='Visit', y=m, connect=True, connect_on=['ptid'], data=plotdf, order=['D0', 'D3', 'D56', 'D59', 'D63', 'D112', 'D168'])
            plt.title(f'{trt}')
            pdf.savefig(figh)