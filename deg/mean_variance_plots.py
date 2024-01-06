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
Average longitudunal plots for modules.
"""

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages


#sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 12

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
# modules_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned_longform.csv')
modules_fn = opj(project_folder, 'wgcna', 'modules_longform_weights.csv')
out_folder = opj(project_folder, 'module_results')
cts_fn = opj(project_folder, 'log_normalized_counts.csv')
res_fn = opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_res.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

cts_df = pd.read_csv(cts_fn)
rx_df = pd.read_csv(rx_fn)
modules_df = pd.read_csv(modules_fn)
res_df = pd.read_csv(res_fn)

sig_df = res_df.loc[(res_df['variable'] == 'day') & (res_df['cohort'] == '3,4') & res_df['sig']]
sig_df = pd.merge(sig_df, modules_df, how='left', on='gene')
sig_df = sig_df.dropna(subset=['module'])

degs = sig_df['gene'].unique().tolist()
all_genes = res_df['gene'].unique().tolist()

cts = cts_df.set_index('Unnamed: 0')

samples = pd.DataFrame(cts.columns, columns=['sampleid'])
samples = samples.assign(ptid=samples['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                         day=samples['sampleid'].map(lambda s: s.split('_')[-1]))
samples = pd.merge(samples, rx_df, how='left', on='ptid')

samps34 = samples.loc[samples['Treatment_Group'].isin(trt34)]

ind_cols = ['sampleid', 'ptid', 'day', 'Treatment_Group', 'Treatment', 'pubid', 'sex']
cts = pd.merge(samps34, cts_df.set_index('Unnamed: 0').T, left_on='sampleid', right_index=True, how='left')

deg_cts = cts.set_index(ind_cols)[degs].stack().reset_index().rename({'level_7':'gene', 0:'GEX'}, axis=1)
deg_cts = pd.merge(deg_cts, modules_df, how='left', on='gene')

all_cts = cts.set_index(ind_cols)[all_genes].stack().reset_index().rename({'level_7':'gene', 0:'GEX'}, axis=1)

tmp = all_cts.groupby('gene')['GEX'].agg([np.mean, np.std, np.var, lambda v: np.sqrt(np.std(v))])
tmp_deg = deg_cts.groupby(['module', 'gene', 'weight'])['GEX'].agg([np.mean, np.var, np.std, lambda v: np.sqrt(np.std(v))])

with PngPdfPages(opj(out_folder, f'mean_variance_plots.pdf')) as pdf:
    for mod in modules_df['module'].unique():
        dat = tmp_deg.loc[mod].reset_index()
        dat = dat.assign(weight=dat['weight'] / dat['weight'].max())

        figh = plt.figure()
        sns.regplot(x='mean', y='<lambda_0>', data=tmp.sample(frac=0.2, random_state=1), lowess=True,
                    line_kws={'color':'black'},
                    scatter_kws={'alpha':0.2, 'color':'gray', 'zorder':-10})
        sns.regplot(x='mean', y='<lambda_0>', data=dat, lowess=True, line_kws={'color':mod}, scatter_kws={'alpha':0.3, 'color':mod, 'zorder':1})
        plt.ylabel('sqrt(standard deviation)')
        plt.xlabel('mean(log-cpm)')
        plt.ylim((0.25, 2.5))
        plt.xlim((-3, 15))
        plt.title(f'{mod.title()} DEGs (n={dat.shape[0]})')
        plt.legend([plt.Circle(1, color=mod), plt.Circle(1, color='gray')], [f'{mod} DEGs (n={dat.shape[0]})', f'All genes (post-filtering, n={tmp.shape[0]})'])
        pdf.savefig(figh)
        plt.close(figh)

        figh = plt.figure()
        sns.regplot(x='mean', y='<lambda_0>', data=tmp.sample(frac=0.2, random_state=1), lowess=True, line_kws={'color':'black'}, scatter_kws={'alpha':0.2, 'color':'gray', 'zorder':-10})
        sns.scatterplot(x='mean', y='<lambda_0>', hue='weight', data=dat)
        plt.ylabel('sqrt(standard deviation)')
        plt.xlabel('mean(log-cpm)')
        plt.ylim((0.25, 2.5))
        plt.xlim((-3, 15))
        plt.title(f'{mod.title()} DEGs (n={dat.shape[0]})')
        pdf.savefig(figh)
        plt.close(figh)

        figh = plt.figure()
        sns.regplot(x='mean', y='weight', data=dat, lowess=True, line_kws={'color':mod}, scatter_kws={'alpha':1, 'color':mod, 'zorder':-10})
        plt.ylabel('eigengene weight (relative)')
        plt.xlabel('mean(log-cpm)')
        plt.xlim((-3, 15))
        plt.ylim((0, 1))
        plt.title(f'{mod.title()} DEGs (n={dat.shape[0]})')
        pdf.savefig(figh)
        plt.close(figh)