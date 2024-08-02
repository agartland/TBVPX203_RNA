import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#import seaborn as sns
import itertools
from functools import partial
import sys
import os
from scipy import stats
import seaborn as sns

from os.path import join as opj

from fg_shared import *

"""SUMMARY:
Heatmap for first figure of manuscript
"""

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from myboxplot import swarmbox
from cluster_alignment import align_clusters

#sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 12

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned_longform.csv')
out_folder = opj(project_folder, 'module_results')
cts_fn = opj(project_folder, 'log_normalized_counts.csv')
#res_fn = opj(project_folder, 'agregated_results_2023-MAR-15.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')
gsea_fn = opj(project_folder, 'gsea', 'agregated_gsea_2024-MAR-14.csv')

#res_df = pd.read_csv(res_fn)
cts_df = pd.read_csv(cts_fn)
rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)

treatments = ['2 µg ID93 + 2 µg GLA-SE',
              '10 µg ID93 + 2 µg GLA-SE',
              '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
              '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)',
              'Placebo']

trt_colors = {'2 µg ID93 + 2 µg GLA-SE':'#00703c',
              '2 µg ID93 + 2 µg GLA-SE (2-dose)':'#00703c',
              '10 µg ID93 + 2 µg GLA-SE':'#00549f',
              '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)':'#ee2e24',
              '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)':'#936fb1',
              '2 µg ID93 + 5 µg GLA-SE (3-dose)':'#ee2e24',
              '2 µg ID93 + 5 µg GLA-SE (2-dose)':'#936fb1',
              'Placebo':'#009fc3'}

# trt_colors = {t:c for t,c in zip(treatments, mpl.cm.tab10.colors)}

trt34 = ['2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)']

cts = cts_df.set_index('Unnamed: 0')
samples = pd.DataFrame(cts.columns, columns=['sampleid'])
samples = samples.assign(ptid=samples['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                         day=samples['sampleid'].map(lambda s: s.split('_')[-1]))
samples = pd.merge(samples, rx_df, how='left', on='ptid')

day_ticks = [0, 1, 2, 3, 4, 5, 6]
day_labels = ['0', '3', '56', '59', '63', '112', '168']
day_map = {k:v for k,v in zip(day_labels, day_ticks)}

mod_groups = [('blue', 'black', 'brown', 'red','green', 'turquoise', 'yellow', 'pink'),
              ('blue',  'brown', 'red'),
              ('green', ),
              ('yellow','pink', 'black','turquoise')]

modnames = {'turquoise':'COREG', 'brown':'BCELL', 'blue':'MITOSIS',
            'pink':'EOS', 'yellow':'IFN-I', 'green':'NEUTRO-I',
            'red':'NEUTRO-II', 'black':'MONO'}
modcolors = {v:k for k,v in modnames.items()}

mod_order = ['IFN-I','EOS','MONO', 'NEUTRO-I',
              'MITOSIS','BCELL', 'NEUTRO-II', 'COREG']
modules = mod_groups[0]

avg = []
for mod in modules:
    genes = modules_df.loc[modules_df['module'] == mod, 'gene'].tolist()
    tmp = cts.loc[genes].mean(axis=0).T.reset_index()
    tmp.columns = ['sampleid', 'score']
    tmp = tmp.assign(module=mod)
    avg.append(tmp)
avg = pd.concat(avg, axis=0)


mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize'] = 9
mpl.rcParams['ytick.labelsize'] = 9
mpl.rcParams['axes.labelsize'] = 10

gene_order = modules_df.sort_values(by='module')['gene'].values
plot_df = pd.merge(samples, cts.loc[gene_order].T, how='left', left_on='sampleid', right_index=True)
# plot_df = plot_df.loc[plot_df['day'].isin(['0', '3', '56', '59', '63']) & plot_df['Treatment_Group'].isin(trt34 + ['Placebo'])]
# plot_df = plot_df.groupby(['Treatment_Group', 'day'])[gene_order].agg('mean').T
delta_df = plot_df.set_index(['pubid', 'day', 'Treatment_Group'])[gene_order].stack().reset_index().rename({'level_3':'Gene', 0:'ct'}, axis=1)

plot_df = plot_df.loc[plot_df['day'].isin(['0', '3', '56', '59', '63']) & plot_df['Treatment_Group'].isin(trt34)]
plot_df.loc[:, gene_order] = plot_df.loc[:, gene_order].apply(stats.zscore, axis=0)
plot_df = plot_df.groupby(['day'])[gene_order].agg('mean').T


delta_df = delta_df.assign(zct=delta_df.groupby('Gene')['ct'].transform(partial(stats.mstats.zscore, nan_policy='omit')))
tmp0 = delta_df.loc[delta_df['day'] == '0']
tmp56 = delta_df.loc[delta_df['day'] == '56']
delta_df = pd.merge(delta_df, tmp0[['pubid', 'Gene', 'zct']], how='left', on=['Gene', 'pubid'], suffixes=('', '_0'))
delta_df = pd.merge(delta_df, tmp56[['pubid', 'Gene', 'zct']], how='left', on=['Gene', 'pubid'], suffixes=('', '_56'))

delta_df = delta_df.assign(delta_0=delta_df['zct'] - delta_df['zct_0'],
                         delta_56=delta_df['zct'] - delta_df['zct_56'])

delta_df = delta_df.loc[delta_df['day'].isin(['3', '59', '63']) & delta_df['Treatment_Group'].isin(trt34)]
delta_56 = delta_df.groupby(['day', 'Gene'])['delta_56'].agg('mean').unstack('Gene').T
delta_0 = delta_df.groupby(['day', 'Gene'])['delta_0'].agg('mean').unstack('Gene').T
delta_56.loc[:, '3'] = delta_0.loc[:, '3']


def _reduced_labels(ind, keep_ind, min_diff=3):
    """Iteratively remove numbers from the sequence until
    the diffs are all at least min_diff"""
    ind = np.asarray(ind)
    max_rem = 10
    i = 0
    while np.any(np.diff(ind) < min_diff): # and i < max_rem:
        i += 1
        # rem_i = np.argmin(np.diff(ind))
        rem_ind = np.argsort(np.diff(ind))
        rem_i = [i for i in rem_ind if not ind[i + 1] in keep_ind][0]
        if np.diff(ind)[rem_i] >= min_diff:
            """This means you have to remove the keep ind"""
            rem_i = np.argmin(np.diff(ind))
            # print(f'Removing keeper: {ind[rem_i + 1]}')
        # print(ind[rem_i + 1], gene_order[ind[rem_i + 1]])
        ind = np.concatenate((ind[:rem_i+1], ind[rem_i+2:]))
    return ind

gsea_res = pd.read_csv(gsea_fn)
gsea_genes = '/'.join(gsea_res['geneID']).split('/')
callouts = "LILRB5,IL8,CD69,DEFA1,DEFA3,DEFA4,LTF,PRTN3,AOC3,OTX1,SDC2,TP35INP2,MRAS,CDKN1C,LYPD2,LYNX1,VMO1,SFTPD,OAS1,BATF2,GPB5,SERPING1".split(',')
gsea_genes = list(set(gsea_genes + callouts))
callout_ind = [i for i,g in enumerate(gene_order) if g in callouts]

# gene_order = modules_df.sample(frac=1).sort_values(by='module')['gene']
ind = [i for i, g in enumerate(gene_order) if g in gsea_genes]
new_ind = _reduced_labels(ind, keep_ind=callout_ind, min_diff=5)
# print(len(new_ind))


with PngPdfPages(opj(out_folder, f'deg_heatmap.pdf'), dpi=300) as pdf:
    cobj = sns.clustermap(data=plot_df,
                          row_cluster=False,
                          col_cluster=False,
                          row_colors=modules_df.sort_values(by='module')['module'].values,
                          figsize=(9, 15),
                          yticklabels=[],
                          vmin=-0.5,
                          vmax=1)
    plt.sca(cobj.ax_heatmap)
    plt.yticks(new_ind, [gene_order[i] for i in new_ind])
    plt.tight_layout()
    plt.xlabel('Study Day')
    plt.ylabel('Gene Module')
    plt.title('\n'.join(trt34))
    pdf.savefig(cobj.fig)

    """Baseline subtracted plots"""
    cobj = sns.clustermap(data=delta_56.loc[gene_order],
                          row_cluster=False,
                          col_cluster=False,
                          row_colors=modules_df.set_index('gene').loc[gene_order, 'module'].values,
                          figsize=(6, 15),
                          yticklabels=[],
                          cmap=mpl.cm.PuOr_r,
                          vmin=-1.5,
                          vmax=1.5)
    plt.sca(cobj.ax_heatmap)
    plt.yticks(new_ind, [gene_order[i] for i in new_ind])
    plt.tight_layout()
    # plt.yticks(ind, gsea_genes)
    plt.xlabel('Study Day')
    plt.ylabel('')
    # plt.title('\n'.join(trt34))
    pdf.savefig(cobj.fig)