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

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from myboxplot import swarmbox
from adjustwithin import adjustwithin
from adjustwithin import adjustnonnan

#sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 12


project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'modules_longform_weights.csv')
out_folder = opj(project_folder, 'correlation_results')
scores_fn = opj(project_folder, 'wgcna', 'module_norm_cts.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

adata_folder = opj(_fg_data, 'SCRI/TBVPX-203/data/immune_adata')
elisa_fn = opj(adata_folder, 'elisa_analysis2.csv')
ics_fn = opj(adata_folder, 'ics_analysis_23.csv')

rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)
mods = modules_df['module'].unique()

treatments = ['2 µg ID93 + 2 µg GLA-SE',
              '10 µg ID93 + 2 µg GLA-SE',
              '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
              '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)',
              'Placebo']

trt_colors = {'2 µg ID93 + 2 µg GLA-SE':'#00703c',
              '10 µg ID93 + 2 µg GLA-SE':'#00549f',
              '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)':'#ee2e24',
              '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)':'#936fb1',
              'Placebo':'#009fc3'}

# trt_colors = {t:c for t,c in zip(treatments, mpl.cm.tab10.colors)}

trt34 = ['2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)']

comparison_map = {'0_3':'Day 3 vs. 0',
                  '56_59':'Day 59 vs. 56',
                  '56_63':'Day 63 vs. 56'}


"""Load and prep ELISA data"""
elisa = pd.read_csv(elisa_fn)
elisa = elisa.assign(ptid=elisa['ptid'].map(lambda v: f'{v:1.0f}'[:3] + '_' + f'{v:1.0f}'[3:]))

keep_days = ['Day 70', 'Day 224']
elisa = elisa.loc[elisa['visitname'].isin(keep_days) & elisa['rx_code'].isin(['T3', 'T4']) & (elisa['analyte'] == 'Total IgG')]

ind_cols = ['ptid', 'visitname']
val_col = 'MEPT'

"""Checks look good"""
# elisa.set_index(ind_cols)[val_col].unstack('visitname')
# pd.merge(elisa, rx_df[['Treatment_Group', 'pubid', 'ptid']], how='left', on='ptid')

elisa = elisa[ind_cols + [val_col]]

"""Load and prep ICS data"""
ics = pd.read_csv(ics_fn, encoding='utf-8', encoding_errors='replace')
ics = ics.assign(ptid=ics['ptid'].map(lambda v: f'{v:1.0f}'[:3] + '_' + f'{v:1.0f}'[3:]),
                 visitname=ics['visitno'].map({70:'Day 70', 224:'Day 224'}))

keep_days = ['Day 70', 'Day 224']
ics = ics.loc[ics['visitname'].isin(keep_days) & ics['Treatment'].isin([3, 4]) & ics['antigen'].isin(['ID93', 'MTB 300'])]

ind_cols = ['ptid', 'visitname', 'tcellsub', 'antigen']
val_col = 'pctpos_adj'

"""Checks look good"""
# ics.set_index(ind_cols)[val_col].unstack(['visitname', 'antigen'])
# pd.merge(ics, rx_df[['Treatment_Group', 'pubid', 'ptid']], how='left', on='ptid')

ics = ics[ind_cols + [val_col]]

"""Load and prep module scores and compute difference variables"""
scores = pd.read_csv(scores_fn)
scores = scores.loc[scores['Treatment'].isin([3, 4])]

ind_cols = ['ptid', 'module', 'day']
val_col = 'eigengene'

scores_0 = scores.loc[scores['day'] == 0]
scores_3 = scores.loc[scores['day'] == 3]
scores_56 = scores.loc[scores['day'] == 56]
scores_5963 = scores.loc[scores['day'].isin([59, 63])]



delta_0 = pd.merge(scores_3[ind_cols + [val_col]],
                scores_0[ind_cols[:2] + [val_col]],
                how='left',
                on=ind_cols[:2],
                suffixes=('', '_0'))
delta_0 = delta_0.assign(delta=delta_0['eigengene'] - delta_0['eigengene_0'],
                         comparison='Day 3 vs. 0')

delta_56 = pd.merge(scores_5963[ind_cols + [val_col]],
                    scores_56[ind_cols[:2] + [val_col]],
                    how='left',
                    on=ind_cols[:2],
                    suffixes=('', '_56'))
delta_56 = delta_56.assign(delta=delta_56['eigengene'] - delta_56['eigengene_56'],
                           comparison=delta_56['day'].map(lambda d: f'Day {d:1.0f} vs. 56'))

cols = ind_cols[:2] + ['delta', 'comparison']
delta = pd.concat((delta_0[cols],
                   delta_56[cols],
                   scores_0.assign(comparison='Day 0', delta=scores_0['eigengene'])[cols],
                   scores_56.assign(comparison='Day 56', delta=scores_56['eigengene'])[cols]), axis=0)

"""Compute pairwise correlations within index-defined groups for ELISA"""
res = []
for elisa_i, elisa_gby in elisa.groupby('visitname'):
    for delta_i, delta_gby in delta.groupby(['comparison', 'module']):
        tmp = pd.merge(elisa_gby, delta_gby, how='inner', on='ptid', validate='1:1')
        tmp_drop = tmp.dropna()
        rho, pvalue = stats.spearmanr(tmp_drop['MEPT'], tmp_drop['delta'])

        res.append(dict(ELISA_day=elisa_i,
                        ELISA_n=elisa_gby.dropna().shape[0],
                        Analyte='Total IgG',
                        ELISA_unit='MEPT',
                        Module=delta_i[1],
                        Comparison=delta_i[0],
                        Module_n=delta_gby.dropna().shape[0],
                        n=tmp_drop.shape[0],
                        rho=rho,
                        pvalue=pvalue))
elisa_res = pd.DataFrame(res)

"""Repeat for ICS"""
res = []
for ics_i, ics_gby in ics.groupby(['visitname', 'tcellsub', 'antigen']):
    for delta_i, delta_gby in delta.groupby(['comparison', 'module']):
        tmp = pd.merge(ics_gby, delta_gby, how='inner', on='ptid', validate='1:1')
        tmp_drop = tmp.dropna()
        rho, pvalue = stats.spearmanr(tmp_drop['pctpos_adj'], tmp_drop['delta'])
        res.append(dict(ICS_day=ics_i[0],
                        tcellsub=ics_i[1],
                        antigen=ics_i[2],
                        ICS_n=ics_gby.dropna().shape[0],
                        ICS_unit='pctps_adj',
                        Module=delta_i[1],
                        Comparison=delta_i[0],
                        Module_n=delta_gby.dropna().shape[0],
                        n=tmp_drop.shape[0],
                        rho=rho,
                        pvalue=pvalue))
ics_res = pd.DataFrame(res)

"""Compute rank tests for deltas t identify significant modules and comparisons"""
def _run_test(i, gby, left_col):
    tmp = gby[[left_col, 'eigengene', 'delta']].dropna()
    w, pvalue = stats.wilcoxon(tmp[left_col], tmp['eigengene'])
    comparison = i[1].replace('Day ', '').replace(' vs. ', '_').split('_')
    out = dict(module=i[0],
                        comparison=i[1],
                        left=comparison[1],
                        right=comparison[0],
                        left_n=gby[left_col].dropna().shape[0],
                        right_n=gby['eigengene'].dropna().shape[0],
                        n=tmp.shape[0],
                        delta=np.mean(tmp['delta']),
                        wstat=w,
                        pvalue=pvalue)
    return out

mod_res = []
for i, gby in delta_0.groupby(['module', 'comparison']):
    out = _run_test(i, gby, left_col='eigengene_0')
    mod_res.append(out)
for i, gby in delta_56.groupby(['module', 'comparison']):
    out = _run_test(i, gby, left_col='eigengene_56')
    mod_res.append(out)
modr = pd.DataFrame(mod_res)

modr = modr.assign(FDRq=adjustnonnan(modr['pvalue'], method='fdr_bh'),
                    FWERp=adjustnonnan(modr['pvalue'], method='holm'))

"""Identify the significant module:comparisons and their individual days as well"""
sig_ind = (modr['pvalue'] < 0.05) & (modr['FDRq'] < 0.10)

keep_test = []
for i, r in modr.loc[sig_ind].iterrows():
    tmp = r['module'] + '_' + r['comparison']
    if not tmp in keep_test:
        keep_test.append(tmp)
    tmp = r['module'] + '_' + f'Day {r["left"]}'
    if not tmp in keep_test:
        keep_test.append(tmp)
    tmp = r['module'] + '_' + f'Day {r["right"]}'
    if not tmp in keep_test:
        keep_test.append(tmp)

"""Apply multiplicity adjustment based on significant module:comparisons"""
def _apply_adj(res, keep_test, keep_ind):
    res = res.assign(mod_comparison=res.apply(lambda r: r['Module'] + '_' + r['Comparison'], axis=1))
    res = res.assign(test=res['mod_comparison'].isin(keep_test) & keep_ind,
                    FDRq=res['pvalue'].copy(),
                    FWERp=res['pvalue'].copy())
    res.loc[~res['test'], 'FDRq'] = np.nan
    res.loc[~res['test'], 'FWERp'] = np.nan

    res = res.assign(FDRq=adjustnonnan(res['FDRq'], method='fdr_bh'),
                     FWERp=adjustnonnan(res['FWERp'], method='holm'))
    #res = res.assign(FDRq=adjustwithin(res, pCol='FDRq', withinCols=['Module'], method='fdr_bh'),
    #                 FWERp=adjustwithin(res, pCol='FWERp', withinCols=['Module'], method='holm'),)
    return res

elisa_keep_ind = elisa_res['ELISA_day'] == 'Day 70'
elisa_res = _apply_adj(elisa_res, keep_test,
                        keep_ind=elisa_keep_ind)
ics_keep_ind = (ics_res['ICS_day'] == 'Day 70') & (ics_res['antigen'] == 'ID93') & (ics_res['tcellsub'] == 'CD4+')
ics_res = _apply_adj(ics_res, keep_test,
                        keep_ind=ics_keep_ind) 

elisa_res.to_csv(opj(out_folder, 'ELISA_module_correlations.csv'), index=False)
ics_res.to_csv(opj(out_folder, 'ICS_module_correlations.csv'), index=False)
modr.to_csv(opj(out_folder, 'module_signedrank_tests.csv'), index=False)


with PngPdfPages(opj(out_folder, f'correlation_heatmaps.pdf')) as pdf:
    plot_df = ics_res.loc[ics_res['Comparison'].str.contains('vs.') & ics_keep_ind]
    plot_df = plot_df.set_index(['Module', 'Comparison'])['rho'].unstack('Comparison')
    
    figh = plt.figure(figsize=(6, 7))
    cobj = sns.heatmap(data=plot_df,
                          vmin=-0.5,
                          vmax=0.5,
                          cmap=mpl.cm.PuOr_r,
                          annot=True,
                          fmt='1.2f')
    plt.xlabel('Comparison')
    plt.ylabel('Gene Module')
    plt.title('ICS: ' + '\n'.join(trt34))
    pdf.savefig(figh)

    plot_df = elisa_res.loc[ics_res['Comparison'].str.contains('vs.') & elisa_keep_ind]
    plot_df = plot_df.set_index(['Module', 'Comparison'])['rho'].unstack('Comparison')
    figh = plt.figure(figsize=(6, 7))
    cobj = sns.heatmap(data=plot_df,
                          vmin=-0.5,
                          vmax=0.5,
                          cmap=mpl.cm.PuOr_r,
                          annot=True,
                          fmt='1.2f')
    plt.xlabel('Comparison')
    plt.ylabel('Gene Module')
    plt.title('ELISA: ' + '\n'.join(trt34))
    pdf.savefig(figh)
