import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
import seaborn as sns

from os.path import join as opj

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from myboxplot import swarmbox

#sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 8
mpl.rcParams['figure.titlesize'] = 8
mpl.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.titlesize'] = 10

"""
Rank correlations of adaptive immune data and module scores

Computes "delta" correlations for ICS and ELISA using Day 56 as baseline
"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
out_folder = opj(project_folder, 'correlation_results')

elisa_res = pd.read_csv(opj(out_folder, 'ELISA_module_correlations.csv'))
ics_res = pd.read_csv(opj(out_folder, 'ICS_module_correlations.csv'))
wbics_res = pd.read_csv(opj(out_folder, 'WB_ICS_module_correlations.csv'))

trt34 = [     '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
              '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)']

comparison_map = {'0_3':'Day 3 vs. 0',
                  '56_59':'Day 59 vs. 56',
                  '56_63':'Day 63 vs. 56'}

modnames = {'turquoise':'COREG', 'brown':'BCELL', 'blue':'MITOSIS',
            'pink':'EOS', 'yellow':'IFN-I', 'green':'NEUTRO-I',
            'red':'NEUTRO-II', 'black':'MONO'}
modcolors = {v:k for k,v in modnames.items()}

mod_order = ['BCELL', 'MITOSIS', 'NEUTRO-II', 'COREG', 'MONO', 'EOS', 'IFN-I','NEUTRO-I']


# sns.set(font_scale=1)

with PngPdfPages(opj(out_folder, f'correlation_heatmaps_fulltp_ms.pdf')) as pdf:
    keep_day = ['Day 14', 'Day 28', 'Day 70', 'Day 224']
    keep_comp = ['Day 3 vs. 0', 'Day 59 vs. 56', 'Day 63 vs. 56']

    ics_keep_ind = ics_res['Day'].isin(keep_day) & ics_res['Comparison'].isin(keep_comp) & (ics_res['antigen'] == 'ID93') & (ics_res['tcellsub'] == 'CD4+') & ics_res['test']
    wbics_keep_ind = wbics_res['Day'].isin(keep_day) & wbics_res['Comparison'].isin(keep_comp) & (wbics_res['antigen'] == 'ID93') & (wbics_res['tcellsub'] == 'CD4+') & wbics_res['test']
    elisa_keep_ind = elisa_res['Day'].isin(keep_day) & elisa_res['Comparison'].isin(keep_comp) & elisa_res['test']
    
    pics = ics_res.loc[ics_keep_ind].assign(Assay='PBMC-ICS')
    pwbics = wbics_res.loc[wbics_keep_ind].assign(Assay='WB-ICS')
    pelisa = elisa_res.loc[elisa_keep_ind].assign(Assay='ELISA')
    
    cols = ['Assay', 'Day', 'Module', 'Comparison', 'rho', 'part_rho', 'part_pvalue', 'FDRq', 'pvalue']
    plot_df = pd.concat((pics[cols], pwbics[cols], pelisa[cols]), axis=0)

    plot_df = plot_df.assign(Module=plot_df['Module'].map(modnames))
    annot_df = plot_df.assign(sig=['+' if p < 0.05 else '' for p in plot_df['pvalue']])
    
    annot_df = annot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['sig'].unstack(['Assay', 'Day'])

    plot_df = plot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['rho'].unstack(['Assay', 'Day'])
    
    x_order = [ ( 'ELISA',  f'Day 14'),
                (  'ELISA', f'Day 28'),
                ( 'ELISA',  f'Day 70'),
                (  'ELISA', f'Day 224'),
                ('PBMC-ICS',f'Day 14'),
                ('PBMC-ICS',f'Day 28'),
                ('PBMC-ICS',f'Day 70'),
                ('PBMC-ICS',f'Day 224'),
                ( 'WB-ICS', f'Day 14'),
                ( 'WB-ICS', f'Day 28'),
                ( 'WB-ICS', f'Day 70'),
                ( 'WB-ICS', f'Day 224')]

    tmp = plot_df[x_order].loc[mod_order]
    color_bars = [modcolors[m] for m in tmp.index.get_level_values(0)]
    cobj = sns.clustermap(data=tmp.reset_index(level='Module', drop=True),
                          vmin=-0.8,
                          vmax=0.8,
                          cmap=mpl.cm.RdYlGn,
                          annot=annot_df[x_order].loc[mod_order],
                          fmt='s',
                          row_cluster=False,
                          col_cluster=False,
                          row_colors=color_bars,
                          dendrogram_ratio=0.01,
                          colors_ratio=0.13,
                          figsize=(5, 5))
    cobj.ax_heatmap.set_ylabel('')
    cobj.ax_heatmap.set_xlabel('')
    cobj.fig.subplots_adjust(right=0.55)
    cobj.ax_cbar.set_position((0.75, 0.25, 0.03, 0.4))
    cobj.ax_cbar.set_ylabel('Rank correlation')
    pdf.savefig(cobj.figure)


"""Identify the significant module:comparisons and their individual days as well"""
modr = pd.read_csv(opj(out_folder, 'module_signedrank_tests.csv'))

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

def _censor_annot(plot_df):
    test_lambda = lambda i, v: f'{v:1.2f}' if f'{modcolors[i]}_{comp}' in keep_test else ''
    annot = plot_df.copy()
    for comp in ['Day 3 vs. 0', 'Day 59 vs. 56', 'Day 63 vs. 56']:
        annot = annot.assign(**{comp:annot.reset_index().apply(lambda r: test_lambda(r['Module'], r[comp]), axis=1).values})
    return annot


for delta_lab, d_lab in [('', ''), ('_delta', u"\u0394")]:
    with PngPdfPages(opj(out_folder, f'correlation_heatmaps{delta_lab}_D70_D224.pdf')) as pdf:
        keep_comp = [f'Day 224{d_lab}', f'Day 70{d_lab}']
        ics_keep_ind = ics_res['Day'].isin(keep_comp) & (ics_res['antigen'] == 'ID93') & (ics_res['tcellsub'] == 'CD4+')
        wbics_keep_ind = wbics_res['Day'].isin(keep_comp) & (wbics_res['antigen'] == 'ID93') & (wbics_res['tcellsub'] == 'CD4+')
        elisa_keep_ind = elisa_res['Day'].isin(keep_comp)
        base_cols = ['Day 0', 'Day 56']
        # base_cols = []
    
        pics = ics_res.loc[~ics_res['Comparison'].isin(base_cols) & ics_res['test'] & ics_keep_ind].assign(Assay='PBMC-ICS')
        pwbics = wbics_res.loc[~wbics_res['Comparison'].isin(base_cols) & wbics_res['test'] & wbics_keep_ind].assign(Assay='WB-ICS')
        pelisa = elisa_res.loc[~elisa_res['Comparison'].isin(base_cols) & elisa_res['test'] & elisa_keep_ind].assign(Assay='ELISA')
        
        cols = ['Assay', 'Day', 'Module', 'Comparison', 'rho', 'FDRq', 'pvalue']
        plot_df = pd.concat((pics[cols], pwbics[cols], pelisa[cols]), axis=0)

        plot_df = plot_df.assign(Module=plot_df['Module'].map(modnames))
        annot_df = plot_df.assign(sig=['+' if p < 0.05 else '' for p in plot_df['pvalue']])
        
        annot_df = annot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['sig'].unstack(['Assay', 'Day'])

        plot_df = plot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['rho'].unstack(['Assay', 'Day'])
        
        x_order = [ ( 'ELISA',  f'Day 70{d_lab}'),
                    (  'ELISA', f'Day 224{d_lab}'),
                    ('PBMC-ICS',f'Day 70{d_lab}'),
                    ('PBMC-ICS',f'Day 224{d_lab}'),
                    ( 'WB-ICS', f'Day 70{d_lab}'),
                    ( 'WB-ICS', f'Day 224{d_lab}')]
        
        tmp = plot_df[x_order].loc[mod_order]
        color_bars = [modcolors[m] for m in tmp.index.get_level_values(0)]
        cobj = sns.clustermap(data=tmp.reset_index(level='Module', drop=True),
                              vmin=-0.8,
                              vmax=0.8,
                              cmap=mpl.cm.RdYlGn,
                              annot=annot_df[x_order].loc[mod_order],
                              fmt='s',
                              row_cluster=False,
                              col_cluster=False,
                              row_colors=color_bars,
                              dendrogram_ratio=0.01,
                              colors_ratio=0.13,
                              figsize=(3, 6))
        cobj.ax_heatmap.set_ylabel('')
        cobj.ax_heatmap.set_xlabel('')
        cobj.fig.subplots_adjust(right=0.4)
        cobj.ax_cbar.set_position((0.78, 0.25, 0.03, 0.4))
        cobj.ax_cbar.set_ylabel('Rank correlation')
        pdf.savefig(cobj.figure)


for delta_lab, d_lab in [('', ''), ('_delta', u"\u0394")]:
    with PngPdfPages(opj(out_folder, f'correlation_heatmaps{delta_lab}_fulltp.pdf')) as pdf:
        keep_comp = ['Day 0', f'Day 14{d_lab}', f'Day 28{d_lab}', f'Day 224{d_lab}', f'Day 70{d_lab}']
        ics_keep_ind = ics_res['Day'].isin(keep_comp) & (ics_res['antigen'] == 'ID93') & (ics_res['tcellsub'] == 'CD4+')
        wbics_keep_ind = wbics_res['Day'].isin(keep_comp) & (wbics_res['antigen'] == 'ID93') & (wbics_res['tcellsub'] == 'CD4+')
        elisa_keep_ind = elisa_res['Day'].isin(keep_comp)
        base_cols = ['Day 0', 'Day 56']
        # base_cols = []
        
        """test variable excludes Day 14 and 28 because they werent primary so we want to add them but only for the
        primary comparisons eg Day 3 IFN. Here we find the 13 primary comparions (the rows of heatmap) and
        apply them to day 14 and 28 too"""
        pics = ics_res.loc[~ics_res['Comparison'].isin(base_cols) & ics_keep_ind].assign(Assay='PBMC-ICS')
        pwbics = wbics_res.loc[~wbics_res['Comparison'].isin(base_cols) & wbics_keep_ind].assign(Assay='WB-ICS')
        pelisa = elisa_res.loc[~elisa_res['Comparison'].isin(base_cols) & elisa_keep_ind].assign(Assay='ELISA')

        primary13_tests = pics.loc[pics['test'], ['Module', 'Comparison']].drop_duplicates().apply(lambda r: '|'.join(r), axis=1).tolist()
        pics.loc[pics[['Module', 'Comparison']].apply(lambda r: '|'.join(r), axis=1).isin(primary13_tests), 'test'] = True
        pwbics.loc[pwbics[['Module', 'Comparison']].apply(lambda r: '|'.join(r), axis=1).isin(primary13_tests), 'test'] = True
        pelisa.loc[pelisa[['Module', 'Comparison']].apply(lambda r: '|'.join(r), axis=1).isin(primary13_tests), 'test'] = True

        pics = pics.loc[pics['test']]
        pwbics = pwbics.loc[pwbics['test']]
        pelisa = pelisa.loc[pelisa['test']]
        
        cols = ['Assay', 'Day', 'Module', 'Comparison', 'rho', 'FDRq', 'pvalue']
        plot_df = pd.concat((pics[cols], pwbics[cols], pelisa[cols]), axis=0)

        plot_df = plot_df.assign(Module=plot_df['Module'].map(modnames))
        annot_df = plot_df.assign(sig=['+' if p < 0.05 else '' for p in plot_df['pvalue']])
        
        annot_df = annot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['sig'].unstack(['Assay', 'Day'])

        plot_df = plot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['rho'].unstack(['Assay', 'Day'])
        
        x_order = [ ( 'ELISA',  f'Day 0'),
                    ( 'ELISA',  f'Day 14{d_lab}'),
                    ( 'ELISA',  f'Day 28{d_lab}'),
                    ( 'ELISA',  f'Day 70{d_lab}'),
                    ( 'ELISA',  f'Day 224{d_lab}'),
                    ('PBMC-ICS',f'Day 0'),
                    ('PBMC-ICS',f'Day 14{d_lab}'),
                    ('PBMC-ICS',f'Day 28{d_lab}'),
                    ('PBMC-ICS',f'Day 70{d_lab}'),
                    ('PBMC-ICS',f'Day 224{d_lab}'),
                    ( 'WB-ICS', f'Day 0'),
                    ( 'WB-ICS', f'Day 14{d_lab}'),
                    ( 'WB-ICS', f'Day 28{d_lab}'),
                    ( 'WB-ICS', f'Day 70{d_lab}'),
                    ( 'WB-ICS', f'Day 224{d_lab}')]
        
        tmp = plot_df[x_order].loc[mod_order]
        color_bars = [modcolors[m] for m in tmp.index.get_level_values(0)]
        cobj = sns.clustermap(data=tmp.reset_index(level='Module', drop=True),
                              vmin=-0.8,
                              vmax=0.8,
                              cmap=mpl.cm.RdYlGn,
                              annot=annot_df[x_order].loc[mod_order],
                              fmt='s',
                              row_cluster=False,
                              col_cluster=False,
                              row_colors=color_bars,
                              dendrogram_ratio=0.01,
                              colors_ratio=0.13,
                              figsize=(5, 6))
        cobj.ax_heatmap.set_ylabel('')
        cobj.ax_heatmap.set_xlabel('')
        cobj.fig.subplots_adjust(right=0.4)
        cobj.ax_cbar.set_position((0.78, 0.25, 0.03, 0.4))
        cobj.ax_cbar.set_ylabel('Rank correlation')
        pdf.savefig(cobj.figure)

with PngPdfPages(opj(out_folder, f'correlation_heatmaps_day14_28.pdf')) as pdf:
    keep_day = [f'Day 14', f'Day 28']
    keep_comp = [f'Day 0', f'Day 3 vs. 0']

    ics_keep_ind = ics_res['Day'].isin(keep_day) & ics_res['Comparison'].isin(keep_comp) & (ics_res['antigen'] == 'ID93') & (ics_res['tcellsub'] == 'CD4+')
    wbics_keep_ind = wbics_res['Day'].isin(keep_day) & wbics_res['Comparison'].isin(keep_comp) & (wbics_res['antigen'] == 'ID93') & (wbics_res['tcellsub'] == 'CD4+')
    elisa_keep_ind = elisa_res['Day'].isin(keep_day) & elisa_res['Comparison'].isin(keep_comp)
    
    pics = ics_res.loc[ics_keep_ind].assign(Assay='PBMC-ICS')
    pwbics = wbics_res.loc[wbics_keep_ind].assign(Assay='WB-ICS')
    pelisa = elisa_res.loc[elisa_keep_ind].assign(Assay='ELISA')
    
    cols = ['Assay', 'Day', 'Module', 'Comparison', 'rho', 'FDRq', 'pvalue']
    plot_df = pd.concat((pics[cols], pwbics[cols], pelisa[cols]), axis=0)

    plot_df = plot_df.assign(Module=plot_df['Module'].map(modnames))
    annot_df = plot_df.assign(sig=['+' if p < 0.05 else '' for p in plot_df['pvalue']])
    
    annot_df = annot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['sig'].unstack(['Assay', 'Day'])

    plot_df = plot_df.set_index(['Assay', 'Day', 'Module', 'Comparison'])['rho'].unstack(['Assay', 'Day'])
    
    x_order = [ ( 'ELISA',  f'Day 14'),
                (  'ELISA', f'Day 28'),
                ('PBMC-ICS',f'Day 14'),
                ('PBMC-ICS',f'Day 28'),
                ( 'WB-ICS', f'Day 14'),
                ( 'WB-ICS', f'Day 28')]
    mod_order3 = ['MONO', 'EOS', 'IFN-I','NEUTRO-I']
    tmp = plot_df[x_order].loc[mod_order]
    color_bars = [modcolors[m] for m in tmp.index.get_level_values(0)]
    cobj = sns.clustermap(data=tmp.reset_index(level='Module', drop=True),
                          vmin=-0.8,
                          vmax=0.8,
                          cmap=mpl.cm.RdYlGn,
                          annot=annot_df[x_order].loc[mod_order],
                          fmt='s',
                          row_cluster=False,
                          col_cluster=False,
                          row_colors=color_bars,
                          dendrogram_ratio=0.01,
                          colors_ratio=0.13,
                          figsize=(3, 4))
    cobj.ax_heatmap.set_ylabel('')
    cobj.ax_heatmap.set_xlabel('')
    cobj.fig.subplots_adjust(right=0.4)
    cobj.ax_cbar.set_position((0.72, 0.25, 0.03, 0.4))
    cobj.ax_cbar.set_ylabel('Rank correlation')
    pdf.savefig(cobj.figure)


for tp in ['Day 70', 'Day 224']:
    for delta_lab, d_lab in [('', ''), ('_delta', u"\u0394")]:
        ics_keep_ind = (ics_res['Day'] == f'{tp}{d_lab}') & (ics_res['antigen'] == 'ID93') & (ics_res['tcellsub'] == 'CD4+')
        wbics_keep_ind = (wbics_res['Day'] == f'{tp}{d_lab}') & (wbics_res['antigen'] == 'ID93') & (wbics_res['tcellsub'] == 'CD4+')
        elisa_keep_ind = elisa_res['Day'] == f'{tp}{d_lab}'
        with PngPdfPages(opj(out_folder, f'correlation_heatmaps{delta_lab}_{tp}.pdf')) as pdf:
            plot_df = ics_res.loc[ics_res['Comparison'].str.contains('vs.') & ics_keep_ind]
            plot_df = plot_df.assign(Module=plot_df['Module'].map(modnames))
            plot_df = plot_df.set_index(['Module', 'Comparison'])['rho'].unstack('Comparison').loc[mod_order]
            
            figh = plt.figure(figsize=(6, 6))
            axh = figh.add_axes([0.15, 0.15, 0.7, 0.7])
            cobj = sns.heatmap(data=plot_df,
                                  vmin=-0.5,
                                  vmax=0.5,
                                  cmap=mpl.cm.coolwarm,
                                  annot=_censor_annot(plot_df),
                                  fmt='s',
                                  ax=axh)
            plt.xlabel('Comparison')
            plt.ylabel('Gene Module')
            plt.title(f'PBMC-ICS {tp}{d_lab}: ' + '\n'.join(trt34))
            pdf.savefig(figh)

            plot_df = wbics_res.loc[ics_res['Comparison'].str.contains('vs.') & wbics_keep_ind]
            plot_df = plot_df.assign(Module=plot_df['Module'].map(modnames))
            plot_df = plot_df.set_index(['Module', 'Comparison'])['rho'].unstack('Comparison').loc[mod_order]
            
            figh = plt.figure(figsize=(6, 6))
            axh = figh.add_axes([0.15, 0.15, 0.7, 0.7])
            cobj = sns.heatmap(data=plot_df,
                                  vmin=-0.5,
                                  vmax=0.5,
                                  cmap=mpl.cm.coolwarm,
                                  annot=_censor_annot(plot_df),
                                  fmt='s',
                                  ax=axh)
            plt.xlabel('Comparison')
            plt.ylabel('Gene Module')
            plt.title(f'WB ICS {tp}{d_lab}: ' + '\n'.join(trt34))
            pdf.savefig(figh)

            plot_df = elisa_res.loc[ics_res['Comparison'].str.contains('vs.') & elisa_keep_ind]
            plot_df = plot_df.assign(Module=plot_df['Module'].map(modnames))
            plot_df = plot_df.set_index(['Module', 'Comparison'])['rho'].unstack('Comparison').loc[mod_order]
            figh = plt.figure(figsize=(6, 6))
            axh = figh.add_axes([0.15, 0.15, 0.7, 0.7])
            cobj = sns.heatmap(data=plot_df,
                                  vmin=-0.5,
                                  vmax=0.5,
                                  cmap=mpl.cm.coolwarm,
                                  annot=_censor_annot(plot_df),
                                  fmt='s',
                                  ax=axh)
            plt.xlabel('Comparison')
            plt.ylabel('Gene Module')
            plt.title(f'ELISA {tp}{d_lab}: ' + '\n'.join(trt34))
            pdf.savefig(figh)

