import pandas as pd
import numpy as np

import sys
import os

from os.path import join as opj
from glob import glob

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from adjustwithin import adjustwithin

"""SUMMARY:
Create an aggregated file of the regression results.
It includes all results bu with indicators for direction and significance.
"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')

"""Concatenate all the results files for the model ~day + sex + (1|ptid)
NOTE: I confirmed that the multiplicity adjustment is within variable and comparison, and across genes."""
res = []
for fn in glob(opj(project_folder, 'degs_day_sex', '*_lme.csv')):
    if not 'all4' in fn:
        tmp = pd.read_csv(fn)
        fn_bits = fn.split(os.path.sep)[-1].split('_')
        tmp = tmp.assign(comparison=f'{fn_bits[1]}_{fn_bits[2]}',
                         Comparison=f'Day {fn_bits[2]} vs. {fn_bits[1]}',
                         cohort='3,4',
                         model='~day + sex')
        res.append(tmp)
    else:
        tmp = pd.read_csv(fn)
        tmp = tmp.assign(comparison=f'0_3',
                         Comparison=f'Day 3 vs. 0',
                         cohort='all4',
                         model='~day + sex')
        if 'adj' in fn:
            tmp = tmp.assign(model='~day + sex + adjuvant')
        res.append(tmp)

res_df = pd.concat(res, axis=0)
# res_df = res_df.loc[res_df['variable'] == 'day']
res_df = res_df.rename({'pval':'pvalue', 'estimate':'log2FC'}, axis=1)

sig_ind = (res_df['pvalue'] < 0.05) & (res_df['FDR'] < 0.2) & (np.abs(res_df['log2FC']) > 0.5)

res_df = res_df.assign(logP=res_df['pvalue'].map(np.log10),
                        Direction=res_df['log2FC'].map(lambda v: {True:'UP', False:'DOWN'}[v>0]),
                        sig=sig_ind)
degs = res_df.loc[res_df['sig'] & (res_df['variable'] == 'day') & (res_df['model'] == '~day + sex') & (res_df['cohort'] == '3,4')]
res_df.to_csv(opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_res.csv'))


"""Aggregate sex:day interaction models, testing only genes that were DEGs"""
res = []
for fn in glob(opj(project_folder, 'degs_day_sex_interaction', '*_lme.csv')):
    if not 'all4' in fn:
        tmp = pd.read_csv(fn)
        fn_bits = fn.split(os.path.sep)[-1].split('_')
        tmp = tmp.assign(comparison=f'{fn_bits[1]}_{fn_bits[2]}',
                         Comparison=f'Day {fn_bits[2]} vs. {fn_bits[1]}',
                         cohort='3,4',
                         model='~day + sex + day:sex')
        res.append(tmp)
    else:
        tmp = pd.read_csv(fn)
        tmp = tmp.assign(comparison=f'0_3',
                         Comparison=f'Day 3 vs. 0',
                         cohort='all4',
                         model='~day + sex + day:sex')
        if 'adj' in fn:
            tmp = tmp.assign(model='~day + sex + day:sex + adjuvant')
        res.append(tmp)

res_df = pd.concat(res, axis=0)
# res_df = res_df.loc[res_df['variable'] == 'day']
res_df = res_df.rename({'pval':'pvalue', 'estimate':'log2FC'}, axis=1)

tmpdf = pd.merge(res_df, degs[['gene', 'comparison', 'pvalue', 'FDR', 'log2FC']], how='left', on=['gene', 'comparison'], suffixes=('', '_deg'))

sig_df = tmpdf.dropna()

sig_df = sig_df.assign(FDR2=adjustwithin(sig_df, pCol='pvalue', withinCols=['model', 'variable', 'comparison'], method='fdr_bh'))
tmp = sig_df.loc[sig_df['variable'] == 'day:sex'].sort_values(by='pvalue')

sig_ind = (res_df['pvalue'] < 0.05) & (res_df['FDR'] < 0.2)

res_df = res_df.assign(logP=res_df['pvalue'].map(np.log10),
                            Direction=res_df['log2FC'].map(lambda v: {True:'UP', False:'DOWN'}[v>0]),
                            sig=sig_ind)
res_df.to_csv(opj(project_folder, 'degs_day_sex_interaction', 'agg_day_sex_interaction_lme_res.csv'))
sig_df.to_csv(opj(project_folder, 'degs_day_sex_interaction', 'agg_deg_day_sex_interaction_lme_res.csv'))

"""Aggregate trt:day interaction models, testing only genes that were DEGs"""
res = []
for fn in glob(opj(project_folder, 'degs_day_sex_trt_interaction', '*_lme.csv')):
    if not 'all4' in fn:
        tmp = pd.read_csv(fn)
        fn_bits = fn.split(os.path.sep)[-1].split('_')
        tmp = tmp.assign(comparison=f'{fn_bits[1]}_{fn_bits[2]}',
                         Comparison=f'Day {fn_bits[2]} vs. {fn_bits[1]}',
                         cohort='3,4',
                         model='~day + sex + trt + trt:day')
        res.append(tmp)
    else:
        tmp = pd.read_csv(fn)
        tmp = tmp.assign(comparison=f'0_3',
                         Comparison=f'Day 3 vs. 0',
                         cohort='all4',
                         model='~day + sex + trt + trt:day')
        if 'adj' in fn:
            tmp = tmp.assign(model='~day + sex + trt + trt:day + adjuvant')
        res.append(tmp)

res_df = pd.concat(res, axis=0)
# res_df = res_df.loc[res_df['variable'] == 'day']
res_df = res_df.rename({'pval':'pvalue', 'estimate':'log2FC'}, axis=1)

tmpdf = pd.merge(res_df, degs[['gene', 'comparison', 'pvalue', 'FDR', 'log2FC']], how='left', on=['gene', 'comparison'], suffixes=('', '_deg'))

sig_df = tmpdf.dropna()

sig_df = sig_df.assign(FDR2=adjustwithin(sig_df, pCol='pvalue', withinCols=['model', 'variable', 'comparison'], method='fdr_bh'))
tmp = sig_df.loc[sig_df['variable'] == 'day:sex'].sort_values(by='pvalue')

sig_ind = (res_df['pvalue'] < 0.05) & (res_df['FDR'] < 0.2)

res_df = res_df.assign(logP=res_df['pvalue'].map(np.log10),
                            Direction=res_df['log2FC'].map(lambda v: {True:'UP', False:'DOWN'}[v>0]),
                            sig=sig_ind)

res_df.to_csv(opj(project_folder, 'degs_day_sex_trt_interaction', 'agg_day_sex_trt_lme_res.csv'))
sig_df.to_csv(opj(project_folder, 'degs_day_sex_trt_interaction', 'agg_deg_day_sex_trt_lme_res.csv'))



"""PLOTS
plot_df = res_df.loc[res_df['variable'].isin(['day', 'sex', 'day:sex'])]
sns.histplot(plot_df.set_index(['gene', 'variable', 'comparison'])['FDR'].unstack(['variable', 'comparison']))

cmap = {v:c for v,c in zip(res_df['variable'].unique(), mpl.cm.tab10.colors)}
plot_df = res_df.loc[res_df['variable'].isin(['day', 'sex'])]
plt.scatter(x='logFDR', y='logFDR2', data=plot_df, c=plot_df['variable'].map(cmap), s=1, marker='.')
xl = plt.xlim()
yl = plt.ylim()
mnmx = [min(xl[0], yl[0]), max(xl[1], yl[1])]
plt.plot(mnmx, mnmx, '--', color='gray')
"""