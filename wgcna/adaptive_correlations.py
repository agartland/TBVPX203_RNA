import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools
from functools import partial
import sys
import os
from scipy import stats

from os.path import join as opj

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from adjustwithin import adjustwithin
from adjustwithin import adjustnonnan
from corrplots import partialcorr

from pngpdf import PngPdfPages
from myboxplot import swarmbox

"""
Rank correlations of adaptive immune data and module scores

Computes "delta" correlations for ICS and ELISA using Day 56 as baseline
"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'modules_longform_weights.csv')
out_folder = opj(project_folder, 'correlation_results')
scores_fn = opj(project_folder, 'wgcna', 'module_norm_cts.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

adata_folder = opj(_fg_data, 'SCRI/TBVPX-203/data/immune_adata')
elisa_fn = opj(adata_folder, 'elisa_analysis2.csv')
ics_fn = opj(adata_folder, 'ics_analysis_23.csv')
wb_ics_fn = opj(adata_folder, 'wb_ics_analysis_24.csv')

rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)
mods = modules_df['module'].unique()

treatments = ['2 µg ID93 + 2 µg GLA-SE',
              '10 µg ID93 + 2 µg GLA-SE',
              '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
              '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)',
              'Placebo']


comparison_map = {'0_3':'Day 3 vs. 0',
                  '56_59':'Day 59 vs. 56',
                  '56_63':'Day 63 vs. 56'}

modnames = {'turquoise':'COREG', 'brown':'BCELL', 'blue':'MITOSIS',
            'pink':'EOS', 'yellow':'IFN-I', 'green':'NEUTRO-I',
            'red':'NEUTRO-II', 'black':'MONO'}
modcolors = {v:k for k,v in modnames.items()}


"""Load and prep ELISA data"""
elisa = pd.read_csv(elisa_fn)

ind_cols = ['ptid', 'rx_code', 'analyte', 'visitname']
val_col = 'MEPT'

elisa = elisa.assign(ptid=elisa['ptid'].map(lambda v: f'{v:1.0f}'[:3] + '_' + f'{v:1.0f}'[3:]))
keep_days = ['Day 0', 'Day 14', 'Day 28', 'Day 56', 'Day 70', 'Day 224']
elisa = elisa.loc[elisa['visitname'].isin(keep_days) & (elisa['analyte'] == 'Total IgG')]
elisa = elisa[ind_cols + [val_col]]

elisa_0 = elisa.loc[elisa['visitname'] == 'Day 56']
elisa = pd.merge(elisa, elisa_0[['ptid', 'analyte', 'rx_code', 'MEPT']], how='left', on=['ptid', 'rx_code', 'analyte'], suffixes=('', '_0'))

elisa = elisa.assign(delta=np.log10(elisa['MEPT']) - np.log10(elisa['MEPT_0']))

keep_days = ['Day 0', 'Day 14', 'Day 28', 'Day 70', 'Day 224']
deltas = elisa.loc[elisa['visitname'].isin(keep_days)]
deltas.loc[:, 'MEPT'] = deltas['delta']
deltas.loc[:, 'visitname'] = deltas['visitname'].map(lambda s: s + u"\u0394")
elisa = pd.concat((elisa.loc[elisa['visitname'].isin(keep_days)],
                    deltas), axis=0).drop('delta', axis=1)

"""Checks look good"""
# elisa.set_index(ind_cols)[val_col].unstack('visitname')
# pd.merge(elisa, rx_df[['Treatment_Group', 'pubid', 'ptid']], how='left', on='ptid')

"""Load and prep ICS data"""
ics = pd.read_csv(ics_fn, encoding='utf-8', encoding_errors='replace')

ind_cols = ['ptid', 'visitname', 'Treatment', 'tcellsub', 'antigen']
val_col = 'pctpos_adj'

ics = ics.assign(ptid=ics['ptid'].map(lambda v: f'{v:1.0f}'[:3] + '_' + f'{v:1.0f}'[3:]),
                 visitname=ics['visitno'].map({70:'Day 70', 56:'Day 56', 224:'Day 224', 0:'Day 0', 14:'Day 14', 28:'Day 28'}))

keep_days = ['Day 0', 'Day 14', 'Day 28', 'Day 56', 'Day 70', 'Day 224']
ics = ics.loc[ics['visitname'].isin(keep_days) & ics['antigen'].isin(['ID93', 'MTB 300'])]
ics = ics[ind_cols + [val_col]]

ics_0 = ics.loc[ics['visitname'] == 'Day 56']
ics = pd.merge(ics, ics_0[['ptid', 'antigen', 'tcellsub', 'Treatment', 'pctpos_adj']],
                how='left',
                on=['ptid', 'Treatment', 'antigen', 'tcellsub'], suffixes=('', '_0'))

ics = ics.assign(delta=ics['pctpos_adj'] - ics['pctpos_adj_0'])

keep_days = ['Day 0', 'Day 14', 'Day 28', 'Day 70', 'Day 224']
deltas = ics.loc[ics['visitname'].isin(keep_days)]
deltas.loc[:, 'pctpos_adj'] = deltas['delta']
deltas.loc[:, 'visitname'] = deltas['visitname'].map(lambda s: s + u"\u0394")
ics = pd.concat((ics.loc[ics['visitname'].isin(keep_days)],
                    deltas), axis=0).drop('delta', axis=1)

"""Checks look good"""
# ics.set_index(ind_cols)[val_col].unstack(['visitname', 'antigen'])
# pd.merge(ics, rx_df[['Treatment_Group', 'pubid', 'ptid']], how='left', on='ptid')


"""Load and prep WB ICS data"""
wbics = pd.read_csv(wb_ics_fn, encoding='utf-8', encoding_errors='replace')

ind_cols = ['ptid', 'visitname', 'Treatment', 'tcellsub', 'antigen']
val_col = 'pctpos_adj'

wbics = wbics.assign(ptid=wbics['ptid'].map(lambda v: f'{v:1.0f}'[:3] + '_' + f'{v:1.0f}'[3:]),
                     visitname=wbics['visitno'].map({70:'Day 70', 56:'Day 56', 224:'Day 224', 0:'Day 0', 14:'Day 14', 28:'Day 28'}))

keep_days = ['Day 0', 'Day 14', 'Day 28', 'Day 56', 'Day 70', 'Day 224']
wbics = wbics.loc[wbics['visitname'].isin(keep_days) & wbics['antigen'].isin(['ID93'])]
wbics = wbics[ind_cols + [val_col]]

wbics_0 = wbics.loc[wbics['visitname'] == 'Day 56']
wbics = pd.merge(wbics, wbics_0[['ptid', 'antigen', 'tcellsub', 'Treatment', 'pctpos_adj']],
                 how='left',
                 on=['ptid', 'Treatment', 'antigen', 'tcellsub'], suffixes=('', '_0'))

wbics = wbics.assign(delta=wbics['pctpos_adj'] - wbics['pctpos_adj_0'])

keep_days = ['Day 0', 'Day 14', 'Day 28', 'Day 70', 'Day 224']
deltas = wbics.loc[wbics['visitname'].isin(keep_days)]
deltas.loc[:, 'pctpos_adj'] = deltas['delta']
deltas.loc[:, 'visitname'] = deltas['visitname'].map(lambda s: s + u"\u0394")
wbics = pd.concat((wbics.loc[wbics['visitname'].isin(keep_days)],
                    deltas), axis=0).drop('delta', axis=1)

"""Checks look good"""
# wbics.set_index(ind_cols)[val_col].unstack(['visitname', 'antigen'])
# pd.merge(wbics, rx_df[['Treatment_Group', 'pubid', 'ptid']], how='left', on='ptid')


"""Load and prep module scores and compute difference variables"""
scores = pd.read_csv(scores_fn)

ind_cols = ['ptid', 'Treatment', 'module', 'day']
val_col = 'eigengene'

scores_0 = scores.loc[scores['day'] == 0]
scores_3 = scores.loc[scores['day'] == 3]
scores_56 = scores.loc[scores['day'] == 56]
scores_5963 = scores.loc[scores['day'].isin([59, 63])]

delta_0 = pd.merge(scores_3[ind_cols + [val_col]],
                scores_0[ind_cols[:3] + [val_col]],
                how='left',
                on=ind_cols[:3],
                suffixes=('', '_0'))
delta_0 = delta_0.assign(delta=delta_0['eigengene'] - delta_0['eigengene_0'],
                         comparison='Day 3 vs. 0')

delta_56 = pd.merge(scores_5963[ind_cols + [val_col]],
                    scores_56[ind_cols[:3] + [val_col]],
                    how='left',
                    on=ind_cols[:3],
                    suffixes=('', '_56'))
delta_56 = delta_56.assign(delta=delta_56['eigengene'] - delta_56['eigengene_56'],
                           comparison=delta_56['day'].map(lambda d: f'Day {d:1.0f} vs. 56'))

cols = ind_cols[:3] + ['delta', 'comparison']
delta = pd.concat((delta_0[cols],
                   delta_56[cols],
                   scores_0.assign(comparison='Day 0', delta=scores_0['eigengene'])[cols],
                   scores_56.assign(comparison='Day 56', delta=scores_56['eigengene'])[cols]), axis=0)

elisa = pd.merge(elisa, rx_df[['Treatment', 'sex', 'ptid']], how='left', on='ptid')
ics = pd.merge(ics, rx_df[['sex', 'ptid']], how='left', on='ptid')
wbics = pd.merge(wbics, rx_df[['sex', 'ptid']], how='left', on='ptid')

elisa = elisa.assign(sex_female=elisa['sex'].map({'female':1, 'male':0}))
ics = ics.assign(sex_female=ics['sex'].map({'female':1, 'male':0}))
wbics = wbics.assign(sex_female=wbics['sex'].map({'female':1, 'male':0}))

"""Save the adata sets to CSV for prediction analysis"""
elisa.to_csv(opj(adata_folder, 'elisa_adata.csv'), index=False)
ics.to_csv(opj(adata_folder, 'ics_adata.csv'), index=False)
wbics.to_csv(opj(adata_folder, 'wbics_adata.csv'), index=False)
delta.to_csv(opj(adata_folder, 'module_adata.csv'), index=False)

"""Drop other treatments for correlation analysis"""
# delta = delta.loc[delta['Treatment'].isin([3, 4])]

"""Only drop placebos"""
delta = delta.loc[delta['Treatment'].isin([1, 2, 3, 4])]

"""Compute pairwise correlations within index-defined groups for ELISA"""
res = []
for elisa_i, elisa_gby in elisa.groupby('visitname'):
    for delta_i, delta_gby in delta.groupby(['comparison', 'module']):
        tmp = pd.merge(elisa_gby, delta_gby, how='inner', on='ptid', validate='1:1')
        tmp_drop = tmp.dropna(subset=['MEPT', 'delta'])
        rho, pvalue = stats.spearmanr(tmp_drop['MEPT'], tmp_drop['delta'])
        part_rho, part_pvalue = partialcorr(tmp_drop['MEPT'],
                                  tmp_drop['delta'],
                                  adjust=[tmp_drop['sex_female']],
                                  method='spearman')

        res.append(dict(Day=elisa_i,
                        day_delta=u"\u0394" in elisa_i[0],
                        Assay_n=elisa_gby.dropna().shape[0],
                        Analyte='Total IgG',
                        ELISA_unit='MEPT',
                        Module=delta_i[1],
                        Comparison=delta_i[0],
                        Module_n=delta_gby.dropna().shape[0],
                        n=tmp_drop.shape[0],
                        rho=rho,
                        pvalue=pvalue,
                        part_rho=part_rho,
                        part_pvalue=part_pvalue))
elisa_res = pd.DataFrame(res)

"""Repeat for ICS"""
res = []
for ics_i, ics_gby in ics.groupby(['visitname', 'tcellsub', 'antigen']):
    for delta_i, delta_gby in delta.groupby(['comparison', 'module']):
        tmp = pd.merge(ics_gby, delta_gby, how='inner', on='ptid', validate='1:1')
        tmp_drop = tmp.dropna(subset=['pctpos_adj', 'delta'])
        rho, pvalue = stats.spearmanr(tmp_drop['pctpos_adj'], tmp_drop['delta'])
        part_rho, part_pvalue = partialcorr(tmp_drop['pctpos_adj'],
                                  tmp_drop['delta'],
                                  adjust=[tmp_drop['sex_female']],
                                  method='spearman')

        res.append(dict(Day=ics_i[0],
                        day_delta=u"\u0394" in ics_i[0],
                        tcellsub=ics_i[1],
                        antigen=ics_i[2],
                        ICS_n=ics_gby.dropna().shape[0],
                        ICS_unit='pctps_adj',
                        Module=delta_i[1],
                        Comparison=delta_i[0],
                        Module_n=delta_gby.dropna().shape[0],
                        n=tmp_drop.shape[0],
                        rho=rho,
                        pvalue=pvalue,
                        part_rho=part_rho,
                        part_pvalue=part_pvalue))
ics_res = pd.DataFrame(res)

"""Repeat for WB ICS"""
res = []
for ics_i, ics_gby in wbics.groupby(['visitname', 'tcellsub', 'antigen']):
    for delta_i, delta_gby in delta.groupby(['comparison', 'module']):
        tmp = pd.merge(ics_gby, delta_gby, how='inner', on='ptid', validate='1:1')
        tmp_drop = tmp.dropna(subset=['pctpos_adj', 'delta'])
        rho, pvalue = stats.spearmanr(tmp_drop['pctpos_adj'], tmp_drop['delta'])
        part_rho, part_pvalue = partialcorr(tmp_drop['pctpos_adj'],
                                  tmp_drop['delta'],
                                  adjust=[tmp_drop['sex_female']],
                                  method='spearman')
        res.append(dict(Day=ics_i[0],
                        day_delta=u"\u0394" in ics_i[0],
                        tcellsub=ics_i[1],
                        antigen=ics_i[2],
                        ICS_n=ics_gby.dropna().shape[0],
                        ICS_unit='pctps_adj',
                        Module=delta_i[1],
                        Comparison=delta_i[0],
                        Module_n=delta_gby.dropna().shape[0],
                        n=tmp_drop.shape[0],
                        rho=rho,
                        pvalue=pvalue,
                        part_rho=part_rho,
                        part_pvalue=part_pvalue))
wbics_res = pd.DataFrame(res)

"""Compute rank tests for deltas to identify significant modules and comparisons"""
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
sig_ind = (modr['pvalue'] < 0.05) & (modr['FDRq'] < 0.05)

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

    # res = res.assign(FDRq=adjustnonnan(res['FDRq'], method='fdr_bh'),
    #                  FWERp=adjustnonnan(res['FWERp'], method='holm'))
    res = res.assign(FDRq=adjustwithin(res, pCol='FDRq', withinCols=['Day'], method='fdr_bh'),
                     FWERp=adjustwithin(res, pCol='FWERp', withinCols=['Day'], method='holm'),)
    return res

keep_days = ['Day 14', 'Day 28', 'Day 70', 'Day 70Δ', 'Day 224', 'Day 224Δ']

elisa_keep_ind = elisa_res['Day'].isin(keep_days)
elisa_res = _apply_adj(elisa_res, keep_test,
                        keep_ind=elisa_keep_ind)
ics_keep_ind = (ics_res['antigen'] == 'ID93') & (ics_res['tcellsub'] == 'CD4+') & ics_res['Day'].isin(keep_days)
ics_res = _apply_adj(ics_res, keep_test,
                        keep_ind=ics_keep_ind)

wbics_keep_ind = (wbics_res['antigen'] == 'ID93') & (wbics_res['tcellsub'] == 'CD4+') & wbics_res['Day'].isin(keep_days)
wbics_res = _apply_adj(wbics_res, keep_test,
                        keep_ind=wbics_keep_ind)

elisa_res.to_csv(opj(out_folder, 'ELISA_module_correlations.csv'), index=False)
ics_res.to_csv(opj(out_folder, 'ICS_module_correlations.csv'), index=False)
wbics_res.to_csv(opj(out_folder, 'WB_ICS_module_correlations.csv'), index=False)
modr.to_csv(opj(out_folder, 'module_signedrank_tests.csv'), index=False)


"""ALSO QUICK PLOTS OF SEX DIFFERENCES FOR ICS< WBICS AND ELISA"""
test_visits = ['Day 14', 'Day 28', 'Day 70', 'Day 224']
visits = test_visits #['Day 14', 'Day 28']
res = []
with PngPdfPages(opj(out_folder, f'sex_adaptive_boxplots.pdf'), dpi=200) as pdf:
    figh = plt.figure()
    tmp = ics.loc[ics['Treatment'].isin([3,4]) & (ics['antigen'] == 'ID93') & (ics['tcellsub'] == 'CD4+')]
    swarmbox(x='visitname', y='pctpos_adj', hue='sex', connect=False, data=tmp.loc[tmp['visitname'].isin(visits)])
    plt.title('PBMC-ICS ID93-specific CD4+ T cells')
    plt.ylabel('% cells w/ \u22652 of IFN\u03b3, IL2, TNF\u03b1')
    pdf.savefig(figh)
    for v in test_visits:
        tmp2 = tmp.loc[tmp['visitname'] == v].dropna()
        s, pvalue = stats.ranksums(tmp2.loc[tmp2['sex']=='female', 'pctpos_adj'], tmp2.loc[tmp2['sex']=='male', 'pctpos_adj'])
        res.append(dict(assay='PBMC-ICS',
                        visit=v,
                        pvalue=pvalue))
    

    figh = plt.figure()
    tmp = wbics.loc[wbics['Treatment'].isin([3,4]) & (wbics['antigen'] == 'ID93') & (wbics['tcellsub'] == 'CD4+')]
    swarmbox(x='visitname', y='pctpos_adj', hue='sex', connect=False, data=tmp.loc[tmp['visitname'].isin(visits)])
    plt.title('WB-ICS ID93-specific CD4+ T cells')
    plt.ylabel('% cells w/ \u22652 of IFN\u03b3, IL2, TNF\u03b1')
    pdf.savefig(figh)
    for v in test_visits:
        tmp2 = tmp.loc[tmp['visitname'] == v].dropna()
        s, pvalue = stats.ranksums(tmp2.loc[tmp2['sex']=='female', 'pctpos_adj'], tmp2.loc[tmp2['sex']=='male', 'pctpos_adj'])
        res.append(dict(assay='WB-ICS',
                        visit=v,
                        pvalue=pvalue))
    

    figh = plt.figure()
    tmp = elisa.loc[elisa['rx_code'].isin(['T3', 'T4']) & (elisa['analyte'] == 'Total IgG')]
    swarmbox(x='visitname', y='MEPT', hue='sex', connect=True, connect_on=['ptid'], data=tmp.loc[tmp['visitname'].isin(visits)])
    plt.yscale('log')
    plt.title('ELISA ID93-specific IgG')
    pdf.savefig(figh)
    for v in test_visits:
        tmp2 = tmp.loc[tmp['visitname'] == v].dropna()
        s, pvalue = stats.ranksums(tmp2.loc[tmp2['sex']=='female', 'MEPT'], tmp2.loc[tmp2['sex']=='male', 'MEPT'])
        res.append(dict(assay='ELISA',
                        visit=v,
                        pvalue=pvalue))

res = pd.DataFrame(res)
res.to_csv(opj(out_folder, 'adaptive_sex_assoc.csv'))
    
    
