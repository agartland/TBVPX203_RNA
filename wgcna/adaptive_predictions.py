import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd
import numpy as np
import itertools
from functools import partial
import sys
import os
from scipy import stats
import seaborn as sns

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn import metrics
from scipy import stats
from sklearn.model_selection import StratifiedKFold, GroupKFold, KFold

from os.path import join as opj

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from myboxplot import swarmbox
from roc import smLogisticRegression, roc_auc_ci, plotROC
from bootstrap_roc import bootstrap_roc

sns.set_style('ticks')
mpl.rcParams['font.size'] = 9
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['axes.labelsize'] = 9

"""
ROC analysis for ELISA, ICS and WB ICS data
"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
out_folder = opj(project_folder, 'correlation_results')
adata_folder = opj(_fg_data, 'SCRI/TBVPX-203/data/immune_adata')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

rx_df = pd.read_csv(rx_fn)

elisa = pd.read_csv(opj(adata_folder, 'elisa_adata.csv'))
ics = pd.read_csv(opj(adata_folder, 'ics_adata.csv'))
wbics = pd.read_csv(opj(adata_folder, 'wbics_adata.csv'))
delta = pd.read_csv(opj(adata_folder, 'module_adata.csv'))

visit_map = {'Day 224':'D224', 'Day 70':'D70', 'Day 224Δ':'D224d', 'Day 70Δ':'D70d'}
comp_map = {'Day 3 vs. 0':'D3d', 'Day 59 vs. 56':'D59d', 'Day 63 vs. 56':'D63d', 'Day 0':'D0', 'Day 56':'D56'}

datasets = {}
elisa = elisa.assign(feat=elisa['visitname'].map(lambda s: 'elisa_' + visit_map[s]))
datasets['elisa'] =  elisa.set_index(['ptid', 'feat'])['MEPT'].unstack('feat').reset_index()

ics = ics.loc[(ics['tcellsub'] == 'CD4+') & (ics['antigen']=='ID93')]
ics = ics.assign(feat=ics['visitname'].map(lambda s: 'ics_' + visit_map[s]))
datasets['ics'] = ics.set_index(['ptid', 'feat'])['pctpos_adj'].unstack('feat').reset_index()

wbics = wbics.loc[(wbics['tcellsub'] == 'CD4+') & (wbics['antigen']=='ID93')]
wbics = wbics.assign(feat=wbics['visitname'].map(lambda s: 'wbics_' + visit_map[s]))
datasets['wbics'] = wbics.set_index(['ptid', 'feat'])['pctpos_adj'].unstack('feat').reset_index()

delta = delta.assign(feat=delta.apply(lambda r: f'gex_{r["module"]}_{comp_map[r["comparison"]]}', axis=1))
datasets['gex'] = delta.set_index(['ptid', 'Treatment', 'feat'])['delta'].unstack('feat').reset_index()

feat_df = pd.merge(datasets['gex'], datasets['elisa'], how='outer', on='ptid')
feat_df = pd.merge(feat_df, datasets['ics'], how='outer', on='ptid')
feat_df = pd.merge(feat_df, datasets['wbics'], how='outer', on='ptid')
feat_df = pd.merge(feat_df, rx_df[['ptid', 'sex', 'Treatment_Group']], how='left', on='ptid')

feat_df = feat_df.loc[~(feat_df['Treatment_Group'] == 'Placebo')]
feat_df = feat_df.assign(low_dose=feat_df['Treatment_Group'].isin(['10 µg ID93 + 2 µg GLA-SE', '2 µg ID93 + 2 µg GLA-SE']).astype(int),
                         three_inj=feat_df['Treatment_Group'].isin(['2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)']).astype(int),
                         sex_female=(feat_df['sex'] == 'female').astype(int))

outcomes = ['elisa_D224', 'elisa_D224d', 'elisa_D70',
           'elisa_D70d', 'ics_D224', 'ics_D224d', 'ics_D70', 'ics_D70d',
           'wbics_D224', 'wbics_D224d', 'wbics_D70', 'wbics_D70d']

outcomes = ['elisa_D224', 'elisa_D70',
           'ics_D224', 'ics_D70', 
           'wbics_D224', 'wbics_D70']

for o in outcomes:
    bin_o = (feat_df[o] >= np.nanmedian(feat_df[o])).astype(float)
    bin_o[feat_df[o].isnull()] = np.nan
    feat_df = feat_df.assign(**{f'{o}_cont':feat_df[o].copy(),
                                f'{o}':bin_o})

def _cv_auc(df, x_cols, y_col, model, n_splits=5):
    stacked = []
    kf = KFold(n_splits=n_splits)
    for i, (train_ind, test_ind) in enumerate(kf.split(df[x_cols], df[y_col])):
        # print(i, len(train_ind), len(test_ind))
        res = mod.fit(X=df[x_cols].iloc[train_ind],
                      y=df[y_col].iloc[train_ind])
        # print(f'{np.sum(mod.coef_ > 0)}/{mod.coef_.shape[1]}')
        proba = res.predict_proba(X=df[x_cols].iloc[test_ind])[:, 1]
        stacked.append(pd.DataFrame({y_col:df[y_col].iloc[test_ind],
                                     'y_pred':proba}))
    stacked = pd.concat(stacked, axis=0)

    # fpr, tpr, thresholds = metrics.roc_curve(stacked['surv'], stacked['y_pred'])
    # pr_threshold = thresholds[np.argmax(tpr - fpr)]
    auc = metrics.roc_auc_score(stacked[y_col], stacked['y_pred'])
    # cmat = metrics.confusion_matrix(stacked['surv'], stacked['y_pred'] > pr_threshold)
    # fres = stats.fisher_exact(cmat)
    lcl, ucl, auc = roc_auc_ci(stacked[y_col].values, stacked['y_pred'].values, alpha=0.05)

    # est, ci = bootstrap_roc(stacked['y_pred'], stacked['surv'], thresholds=100)

    stacked = stacked.assign(variable=' + '.join(x_cols),
                             auc=auc, auc_lcl=lcl, auc_ucl=ucl,
                             n_vars=len(x_cols),
                             n_nz_vars=np.sum(mod.coef_ > 0))
    return stacked

var_0 = [  'gex_black_D0', 'gex_black_D56', # 12
           'gex_blue_D56',
           'gex_brown_D56',
           'gex_green_D0', 'gex_green_D56',
           'gex_pink_D0',  'gex_pink_D56',
           'gex_red_D56', 
           'gex_turquoise_D56',
           'gex_yellow_D0', 'gex_yellow_D56']

var_boost = ['gex_black_D59d','gex_black_D63d', # 9
               'gex_blue_D63d',
               'gex_brown_D63d',
               'gex_green_D59d',
               'gex_pink_D59d',           
               'gex_red_D63d', 
               'gex_turquoise_D63d',
               'gex_yellow_D59d'] 

var_prime = ['gex_black_D3d',  # 4
               'gex_green_D3d',
               'gex_pink_D3d',
               'gex_yellow_D3d']

var_ifn = ['gex_yellow_D3d', 'gex_yellow_D59d'] # 2

base_cols = ['low_dose', 'three_inj', 'sex_female']

"""Try univariate models of each variable"""
mod = LogisticRegression(penalty=None, solver='lbfgs')
results = []
for y_col in outcomes:
    for col in var_boost + var_prime + var_0:
        stacked = _cv_auc(feat_df[[col, y_col]].dropna(), [col], y_col, mod)
        stacked = stacked.assign(outcome=y_col) # variable=stacked['variable'].str.replace('low_dose + three_inj + sex_female', ''))
        stacked = stacked.drop([y_col, 'y_pred'], axis=1).drop_duplicates()
        results.append(stacked)
results = pd.concat(results, axis=0)
res_summ = results.sort_values('auc', ascending=False)
res_summ.to_csv(opj(out_folder, 'univariate_roc_results.csv'), index=False)

"""Try combinations of variables (includes univariate)"""
mod = LogisticRegression(penalty='l1', solver='liblinear', C=1)
results = []
varlists = [var_ifn, var_prime, var_boost, var_prime + var_boost]
for y_col in outcomes:
    for cols in varlists:
        stacked = _cv_auc(feat_df[[y_col] + cols].dropna(), cols, y_col, mod)
        stacked = stacked.assign(outcome=y_col, y_obs=stacked[y_col]).drop(y_col, axis=1)
        results.append(stacked)
l1_results = pd.concat(results, axis=0)
l1_results = l1_results.assign(Set=l1_results['n_vars'].map({9:'Boost', 4:'Prime', 2:'IFN-I', 13:'Prime + Boost'}))

l1_summ = l1_results.groupby(['outcome', 'variable', 'Set', 'auc', 'auc_lcl', 'auc_ucl', 'n_vars'])['n_nz_vars'].agg('mean').reset_index()
l1_summ = l1_summ.sort_values('auc', ascending=False)
l1_summ.to_csv(opj(out_folder, 'l1_roc_results.csv'), index=False)

def _label_fun(s):
    assay, tp = s.split('_')
    assay = {'elisa':'ELISA IgG', 'ics':'PBMC-ICS CD4+ T', 'wbics':'WB-ICS CD4+ T'}[assay]
    return f'{assay}-{tp}'


def _auc_forest_plot(res_summ, colors={}):
    x_col ='auc'
    x = res_summ[x_col].values
    y = list(range(res_summ.shape[0]))[::-1]
    xerr = np.concatenate(( (res_summ[x_col].values - res_summ[x_col + '_lcl'].values)[None, :],
                            (res_summ[x_col + '_ucl'].values - res_summ[x_col].values)[None, :] ), axis=0)
    figh = plt.figure(figsize=(4, 2))
    axh = figh.add_axes([0.25, 0.2, 0.65, 0.65])
    # plt.errorbar(x, y, xerr=xerr, fmt='o', color='k')
    for i in range(x.shape[0]):
        color = colors.get(res_summ['Set'].iloc[i], 'k')
        plt.errorbar(x[i], y[i], xerr=xerr[:, i][:, None], fmt='s', color=color, elinewidth=1, markersize=3)
    yl = (-1, x.shape[0])
    plt.plot([0.5, 0.5], yl, '--', color='gray')
    plt.ylim(yl)
    plt.yticks(y, res_summ['Set'])
    # plt.xticks(np.arange(0.4, 1, 0.05))
    plt.xlim((0.2, 1.02))
    axh.xaxis.set_major_locator(MultipleLocator(0.1))
    axh.xaxis.set_major_formatter('{x:1.2f}')
    axh.xaxis.set_minor_locator(MultipleLocator(0.025))
    plt.xlabel(f'CV-AUROC')
    plt.title(f'{_label_fun(y_col)}', fontsize=9)
    return figh

with PngPdfPages(opj(out_folder, f'ROC_plots.pdf')) as pdf:
    for y_col in  ['elisa_D70','elisa_D224', 'ics_D70', 'ics_D224', 'wbics_D70', 'wbics_D224']:
        lines = ['Prime', 'Boost', 'IFN-I', 'Prime + Boost']

        """
        best = 'Prime + Boost'
        ind = (l1_results['Set'] == best) & (l1_results['outcome'] == y_col)
        est, ci = bootstrap_roc(l1_results.loc[ind, 'y_pred'],
                                l1_results.loc[ind, 'y_obs'],
                                thresholds=100)"""
        figh = _auc_forest_plot(l1_summ.loc[l1_summ['outcome'] == y_col])
        pdf.savefig(figh)
        plt.close(figh)

        figh = plt.figure(figsize=(5, 4))
        axh = figh.add_axes([0.1, 0.12, 0.8, 0.60])
        axh.set_aspect('equal')

        used_colors = {}
        for v, c in zip(lines, mpl.cm.tab10.colors):
            used_colors[v] = c
            ind = (l1_results['Set'] == v) & (l1_results['outcome'] == y_col)
            st = l1_results.loc[ind]
            fpr, tpr, thresholds = metrics.roc_curve(st['y_obs'], st['y_pred'])
            lcl, ucl, auc = roc_auc_ci(st['y_obs'].values, st['y_pred'].values, alpha=0.05)
            label = f'{v} ({auc:1.2f} [{lcl:1.2f} - {ucl:1.2f}])'
            plt.plot(fpr, tpr, color=c, lw=1.5, label=label)
        """plt.fill_between(1 - est['Specificity'],
                         ci['Sensitivity'][:, 1],
                         ci['Sensitivity'][:, 0],
                         color=used_colors[best], alpha=0.3, zorder=-4)"""
        plt.plot([0, 1], [0, 1], '--', color='gray')
        plt.xlim([-0.025, 1.025])
        plt.xlim([-0.025, 1.025])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc='lower center', fontsize=9, bbox_to_anchor=(0.5, 1), title=f'{_label_fun(y_col)} (CV-AUC [95% CI])', title_fontsize=9)
        pdf.savefig(figh)
        plt.close(figh)

with PngPdfPages(opj(out_folder, f'scatter_plots.pdf')) as pdf:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 7))
    g = sns.scatterplot(x='gex_yellow_D59d', y='ics_D224_cont', data=feat_df, color='black', ax=ax1)
    g.axes.set_yscale('log')
    ax1.set_xlabel('IFN-I Day 59 vs. Day 56')
    ax1.set_ylabel('PBMC-ICS CD4+ T cell response\n(% cells w/ \u22652 of IFN\u03b3, IL2, TNF\u03b1)')
    
    g = sns.scatterplot(x='gex_yellow_D3d', y='wbics_D224_cont', data=feat_df, color='black', ax=ax2)
    g.axes.set_yscale('log')
    ax2.set_xlabel('IFN-I Day 3 vs. Day 0')
    ax2.set_ylabel('WB-ICS CD4+ T cell response\n(% cells w/ \u22652 of IFN\u03b3, IL2, TNF\u03b1)')

    pdf.savefig(fig)
