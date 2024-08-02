import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import itertools
import scipy.cluster.hierarchy as sch
import sys
import os

from os.path import join as opj

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
adata_folder = opj(_fg_data, 'SCRI/TBVPX-203/data/immune_adata')

modules_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned_longform.csv')
out_folder = opj(project_folder, 'module_results')
out_file = opj(out_folder, f'gene_by_module_heatmaps.pdf')
cts_fn = opj(project_folder, 'log_normalized_counts.csv')
scores_fn = opj(project_folder, 'wgcna', 'module_norm_cts.csv')
res_fn = opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_res.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

# res_df = pd.read_csv(res_fn)
cts_df = pd.read_csv(cts_fn)
rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)
"""Somehow there are three genes that made the filter for DEG testing
but missed the filter for normalized counts matrix: dropping them here but should be added back to normalized cts"""
modules_df = modules_df.loc[~modules_df['gene'].isin(['MYO5C','DDX11L5', 'LOC101928891'])]

treatments = ['2 µg ID93 + 2 µg GLA-SE', 'Placebo',
               '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)',
               '10 µg ID93 + 2 µg GLA-SE']
trt34 = ['2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)']

modnames = {'turquoise':'COREG', 'brown':'BCELL', 'blue':'MITOSIS',
            'pink':'EOS', 'yellow':'IFN-I', 'green':'NEUTRO-I',
            'red':'NEUTRO-II', 'black':'MONO'}
modcolors = {v:k for k,v in modnames.items()}

corr_threshold = 0.6

res_df = pd.read_csv(res_fn)
sig_df = res_df.loc[(res_df['variable'] == 'day') & (res_df['cohort'] == '3,4') & res_df['sig']]
sig_df = pd.merge(sig_df, modules_df, how='left', on='gene')
sig_df = sig_df.dropna(subset=['module'])
sig_genes = sig_df['gene'].unique().tolist()

cts = cts_df.set_index('Unnamed: 0').loc[sig_genes]

samples = pd.DataFrame(cts.columns, columns=['sampleid'])
samples = samples.assign(ptid=samples['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                         day=samples['sampleid'].map(lambda s: s.split('_')[-1]))
samples = pd.merge(samples, rx_df, how='left', on='ptid')

samps34 = samples.loc[samples['Treatment_Group'].isin(trt34)]

"""Mat has genes as features (columns) and observations as rows"""
mat = cts[samps34['sampleid']].T

"""Using method='pearson' here makes it identical to the empirical covariance matrix"""
corr_mat = mat.corr(method='spearman')

corr_mat.index.name = ''
corr_mat.columns.name = ''

mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6

sorted_genes = {}
all_genes = []
module_order = ['brown', 'blue', 'red', 'turquoise', 'black','pink', 'yellow', 'green']
for m in ['brown', 'blue', 'red', 'turquoise', 'black','pink', 'yellow', 'green']:
    mod_genes = modules_df.loc[modules_df['module'] == m, 'gene'].values
    Z = sch.linkage(corr_mat.loc[mod_genes, :].loc[:, mod_genes].T, 'ward')
    sorted_genes[m] = mod_genes[sch.leaves_list(Z)].tolist()
    all_genes.extend(sorted_genes[m])

ordered_modules_df = modules_df.assign(name=modules_df['module'].map(modnames)).set_index('gene').loc[all_genes].reset_index()
ordered_modules_df.to_csv(modules_fn.replace('.csv', '_ORDERED.csv'), index=False)

def genes2colors(glist):
    return [modules_df.loc[modules_df['gene'] == g, 'module'].iloc[0] for g in glist]

params = dict(vmin=-1, vmax=1, cmap='RdYlGn',
                    col_cluster=False,
                    row_cluster=False)

with PngPdfPages(out_file) as pdf:
    tmp = corr_mat.loc[all_genes, :].loc[:, all_genes]
    g = sns.clustermap(tmp, 
                        yticklabels=False,
                        xticklabels=False,
                        col_colors=genes2colors(tmp.columns),
                        row_colors=genes2colors(tmp.index),
                        figsize=(10, 14),
                        mask=np.abs(tmp) < 0.2,
                        **params)
    pdf.savefig(g.figure)

    tmp = corr_mat.loc[sorted_genes['green'], :].loc[:, sorted_genes['blue'] + sorted_genes['yellow']]
    g = sns.clustermap(tmp,
                    figsize=(5, 4),
                    col_colors=genes2colors(tmp.columns),
                    row_colors=genes2colors(tmp.index),
                    mask=np.abs(tmp) < 0.2,
                    yticklabels='auto',
                    xticklabels='auto',
                    **params)
    pdf.savefig(g.figure)

    tmp = corr_mat.loc[sorted_genes['turquoise'], :].loc[:, sorted_genes['blue'] + sorted_genes['brown'] + sorted_genes['red']]
    g = sns.clustermap(tmp,
                    figsize=(5, 7),
                    col_colors=genes2colors(tmp.columns),
                    row_colors=genes2colors(tmp.index),
                    mask=np.abs(tmp) < 0.2,
                    yticklabels='auto',
                    xticklabels='auto',
                    **params)
    pdf.savefig(g.figure)

    tmp = corr_mat.loc[sorted_genes['turquoise'], :].loc[:, sorted_genes['red']]
    g = sns.clustermap(tmp,
                    figsize=(5, 7),
                    col_colors=genes2colors(tmp.columns),
                    row_colors=genes2colors(tmp.index),
                    yticklabels='auto',
                    xticklabels='auto',
                    **params)
    pdf.savefig(g.figure)

    tmp = corr_mat.loc[sorted_genes['turquoise'], :].loc[:, sorted_genes['brown']]
    g = sns.clustermap(tmp,
                    figsize=(5, 7),
                    col_colors=genes2colors(tmp.columns),
                    row_colors=genes2colors(tmp.index),
                    yticklabels='auto',
                    xticklabels='auto',
                    **params)
    pdf.savefig(g.figure)


"""CORRELATION HEATMAPS OF MODULES"""
scores = pd.read_csv(scores_fn)
delta = pd.read_csv(opj(adata_folder, 'module_adata.csv'))

day_colors = {0:'black', 56:'gray', 3:'gold', 59:'orange', 63:'dodgerblue',
              'Day 0':'black', 'Day 56':'gray', 'Day 3 vs. 0':'gold', 'Day 59 vs. 56':'orange', 'Day 63 vs. 56':'dodgerblue'}

mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10

with PngPdfPages(opj(out_folder, f'module_correlation_heatmaps.pdf')) as pdf:
    plotdf = scores.loc[(scores['Treatment_Group'] != 'Placebo') & scores['day'].isin([0, 3, 56, 59, 63])]
    plotdf = plotdf.assign(Module=[modnames[c] for c in plotdf['module']])
    plotdf = plotdf.set_index(['ptid', 'day', 'Module'])['eigengene'].unstack(['Module', 'day']).corr(method='spearman')

    color_bars = pd.DataFrame({'Module':[modcolors[m] for m in plotdf.columns.get_level_values(0)],
                                'Comparison':[day_colors[d] for d in plotdf.columns.get_level_values(1)]}, index=plotdf.index)
    g = sns.clustermap(plotdf,
                   vmin=-0.8, vmax=0.8, cmap='RdYlGn',
                   row_colors=color_bars,
                   col_colors=color_bars,
                   row_cluster=False,
                   col_cluster=False,
                   dendrogram_ratio=0.01,
                   colors_ratio=0.08)
    pdf.savefig(g.figure)

    plotdf = delta.loc[(delta['Treatment'] != 5) & delta['comparison'].isin(['Day 3 vs. 0','Day 59 vs. 56', 'Day 63 vs. 56'])]
    plotdf = plotdf.assign(Module=[modnames[c] for c in plotdf['module']])
    plotdf = plotdf.set_index(['ptid', 'comparison', 'Module'])['delta'].unstack(['Module', 'comparison']).corr(method='spearman')

    color_bars = pd.DataFrame({'Module':[modcolors[m] for m in plotdf.columns.get_level_values(0)],
                                'Comparison':[day_colors[d] for d in plotdf.columns.get_level_values(1)]}, index=plotdf.index)
    g = sns.clustermap(plotdf,
                   vmin=-0.8, vmax=0.8, cmap='RdYlGn',
                   row_colors=color_bars,
                   col_colors=color_bars,
                   row_cluster=True,
                   col_cluster=True,
                   dendrogram_ratio=0.01,
                   colors_ratio=0.08)
    pdf.savefig(g.figure)


"""HEATMAPS OF GENE EXPRESSION BY IMMUNE CELL TYPE"""
imm_db = pd.read_csv(opj(_fg_data, 'immune_genes_db', 'rna_immune_cell.tsv.zip'), compression='zip', sep='\t')

"""Same data but at the sample level"""
# imm_db2 = pd.read_csv(opj(_fg_data, 'immune_genes_db', 'rna_immune_cell_sample.tsv.zip'), compression='zip', sep='\t')

# g20 = tmp.sum(axis=1).sort_values().head(20)
# sns.barplot(x='Gene name', y='nTPM', hue='Immune cell', data=imm_db.loc[imm_db['Gene name'].isin(g20.index)])
# plt.legend(loc='upper left', bbox_to_anchor=(1,1))

# unstacked = imm_db2.loc[imm_db2['Gene name'].isin(g20.index)].groupby(['Immune cell', 'Gene name'])['nTPM'].agg('mean').unstack('Immune cell')

with PngPdfPages(opj(out_folder, f'celltype_heatmaps.pdf')) as pdf:
    for m in ['brown', 'blue', 'red', 'turquoise', 'black','pink', 'yellow', 'green']:
        glist = sorted_genes[m]
        tpm = imm_db.loc[imm_db['Gene name'].isin(glist)].set_index(['Immune cell', 'Gene name'])['nTPM'].unstack('Immune cell')
        tpm = np.log10(tpm)
        tpm.values[tpm.values < 0] = 0
        cobj = sns.clustermap(tpm, cmap='magma_r', figsize=(5, 10), yticklabels=True)
        pdf.savefig(cobj.figure)