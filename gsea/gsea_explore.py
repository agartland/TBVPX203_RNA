import pandas as pd
import numpy as np
import seaborn as sns
import sys
import os

from os.path import join as opj

from gseapy import Msigdb

from fg_shared import *

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned_longform.csv')
cts_fn = opj(project_folder, 'log_normalized_counts.csv')
res_fn = opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_res.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

cts_df = pd.read_csv(cts_fn)
rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)

cts = cts_df.set_index('Unnamed: 0').loc[modules_df['gene'].values]

samples = pd.DataFrame(cts.columns, columns=['sampleid'])
samples = samples.assign(ptid=samples['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                         day=samples['sampleid'].map(lambda s: s.split('_')[-1]))
samples = pd.merge(samples, rx_df, how='left', on='ptid')

samps34 = samples.loc[samples['Treatment_Group'].isin(trt34)]

mat = cts[samps34['sampleid']].T
corr_mat = mat.corr(method='spearman')

genes_ss = modules_df.loc[modules_df['module'] == 'turquoise', 'gene'].tolist()
genes_ss.remove('MYO5C')

corr_mat_ss = corr_mat.loc[:, genes_ss].loc[genes_ss, :]
mpl.rcParams['ytick.labelsize'] = 7
g = sns.clustermap(corr_mat_ss, cmap='Spectral', vmin=-1, vmax=1)
tmp = g.data2d.loc[:'GRAMD1C', :].loc[:, :'GRAMD1C']
g2 = sns.clustermap(tmp, cmap='Spectral', vmin=-1, vmax=1, figsize=(10, 20))


msig = Msigdb()
# mouse hallmark gene sets
gmt = msig.get_gmt(category='mh.all', dbver="2023.1.Mm"

# list msigdb version you wanna query
msig.list_dbver()
# list categories given dbver.
msig.list_category(dbver="2023.1.Hs") # mouse

print(gmt['HALLMARK_WNT_BETA_CATENIN_SIGNALING'])

names = gp.get_library_name()

gene_sets="./data/genes.gmt",
gene_sets={'A':['gene1', 'gene2',...],
           'B':['gene2', 'gene4',...],
           ...}

enr2 = gp.enrich(gene_list="./tests/data/gene_list.txt", # or gene_list=glist
                 gene_sets=["./tests/data/genes.gmt", "unknown", kegg ], # kegg is a dict object
                 background=None, # or "hsapiens_gene_ensembl", or int, or text file, or a list of genes
                 outdir=None,
                 verbose=True)
enr2.results.head()