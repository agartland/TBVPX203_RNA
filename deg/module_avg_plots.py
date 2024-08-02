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

with PngPdfPages(opj(out_folder, f'delta_module_lineplots.pdf')) as pdf:
    """Focus on Day 0 vs 3 (D0 subtracted) and Day 56, 59, 63 (D56 subtracted)
    Each group separately, showing individuals."""

    for mod in modules:
        #mod = 'black'
        plot_df = pd.merge(samples, avg.loc[avg['module'] == mod], how='left', on='sampleid')
        plot_df = plot_df.assign(Day=plot_df['day'].map(day_map))

        tmp0 = plot_df.loc[plot_df['day'] == '0']
        tmp56 = plot_df.loc[plot_df['day'] == '56']
        plot_df = pd.merge(plot_df, tmp0[['pubid', 'module', 'score']], how='left', on=['module', 'pubid'], suffixes=('', '_0'))
        plot_df = pd.merge(plot_df, tmp56[['pubid', 'module', 'score']], how='left', on=['module', 'pubid'], suffixes=('', '_56'))

        plot_df = plot_df.assign(delta_0=plot_df['score'] - plot_df['score_0'],
                                 delta_56=plot_df['score'] - plot_df['score_56'])
        plot_df = plot_df.assign(Treatment_Group=plot_df['Treatment_Group'].str.replace(' Vaccine Injections', '-dose'))

        sex = 'both'
        # for sex in ['both', 'male', 'female']:
        for keep_treatments in [treatments[:1] + treatments[2:], treatments[2:]]:
            width = len(keep_treatments) * 3
            for zero in ['delta_0', 'delta_56']:
                figh = plt.figure(figsize=(width, 4))
                gs = mpl.gridspec.GridSpec(1, len(keep_treatments),
                                           left=0.15,
                                           wspace=0.25,
                                           bottom=0.15)
                axh = []
                for trti, trt in enumerate(keep_treatments):
                    trt = trt.replace(' Vaccine Injections', '-dose')
                    # axh = figh.add_axes([0.12, 0.12, 0.5, 0.8])
                    if trti == 0:
                        axh.append(figh.add_subplot(gs[0, trti]))
                    else:
                        axh.append(figh.add_subplot(gs[0, trti], sharey=axh[0]))
                    
                    if sex == 'both':
                        ind = plot_df['Treatment_Group'] == trt
                        sex_lab = ''
                    else:
                        ind = (plot_df['Treatment_Group'] == trt) & (plot_df['sex'] == sex)
                        sex_lab = sex
                    
                    tmp = plot_df.loc[ind].set_index(['pubid', 'Day'])[zero].unstack('Day')
                    if zero == 'delta_56':
                        tmp.loc[:, [0, 1]] = np.nan
                        tmp = tmp.drop([0, 1], axis=1)
                    
                    plt.plot(tmp.columns, tmp.values.T, '-', color=trt_colors[trt], label=trt, lw=0.5)
                    plt.plot(tmp.columns, tmp.mean(axis=0), '-s', color=trt_colors[trt], label=trt, lw=2)
                    plt.plot(tmp.columns, [0]*tmp.shape[1], '--', color='gray', zorder=-1)
                    if zero == 'delta_56':
                        plt.xticks(ticks=day_ticks[2:], labels=day_labels[2:])
                    else:
                        plt.xticks(ticks=day_ticks, labels=day_labels)
                    if trti == 0:
                        plt.ylabel(f'{sex_lab.title()}\n{modnames[mod]} Module (normalized counts)')
                    plt.xlabel('Study Day')
                    plt.title(trt.replace('SE (', 'SE\n('), size=10)
                    # plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                pdf.savefig(figh)
                plt.close(figh)


plot_df = pd.merge(samples, avg, how='left', on='sampleid')
# plot_df.to_csv(opj(out_folder, 'module_norm_cts_avg.csv'))

with PngPdfPages(opj(out_folder, f'module_sex2_pubid_lineplots_day3.pdf')) as pdf:
    """Prep scores for summary plot with multiple modules. Mean across genes."""
    trt = trt34
    for m in modules:
        ind = plot_df['Treatment_Group'].isin(trt) & (plot_df['module'] == m)
        ind = ind & plot_df['day'].isin(['0', '3'])
        female = plot_df.loc[ind & (plot_df['sex'] == 'female')]
        male = plot_df.loc[ind & (plot_df['sex'] == 'male')]
        figh = plt.figure(figsize=(4, 4))
        axh = figh.add_axes([0.15, 0.2, 0.5, 0.6])
        swarmbox(x='sex', y='score',
                 hue='day',
                 hue_order=['0', '3'],
                 connect=True, connect_on=['pubid'],
                 data=plot_df.loc[ind],
                 order=['female', 'male'],
                 legend_loc=None,
                 swarm_size=3,
                 box_palette=['lightgray'],
                 swarm_color='k')
        plt.title('\n'.join(trt) + '\n')
        # plt.legend(loc='upper right', title='Day')
        plt.xlabel('Participant Sex')
        plt.ylabel(f'{modnames[m]} Module Score')
        pdf.savefig(figh)
        plt.close(figh)

with PngPdfPages(opj(out_folder, f'module_sex2_pubid_lineplots_day59.pdf')) as pdf:
    """Prep scores for summary plot with multiple modules. Mean across genes."""
    trt = trt34
    trt = trt34[1:]
    for m in modules:
        ind = plot_df['Treatment_Group'].isin(trt) & (plot_df['module'] == m)
        ind = ind & plot_df['day'].isin(['56', '59', '63'])
        #female = plot_df.loc[ind & (plot_df['sex'] == 'female')]
        #male = plot_df.loc[ind & (plot_df['sex'] == 'male')]
        figh = plt.figure(figsize=(5, 4))
        axh = figh.add_axes([0.15, 0.2, 0.5, 0.6])
        swarmbox(x='sex', y='score',
                 hue='day',
                 hue_order=['56', '59', '63'],
                 connect=True, connect_on=['pubid'],
                 data=plot_df.loc[ind],
                 order=['female', 'male'],
                 legend_loc=None,
                 swarm_size=3)
        plt.title('\n'.join(trt))
        # plt.legend(loc='upper right', title='Day')
        plt.xlabel('Participant Sex')
        plt.ylabel(f'{modnames[m]} Module Score')
        pdf.savefig(figh)
        plt.close(figh)


with PngPdfPages(opj(out_folder, f'module_pubid_lineplots.pdf')) as pdf:
    """Prep scores for summary plot with multiple modules. Mean across genes."""
    for trt in [trt34] + [[t] for t in treatments]:
        for m in modules:
            ind = plot_df['Treatment_Group'].isin(trt) & (plot_df['module'] == m)
            figh = plt.figure(figsize=(6, 4))
            axh = figh.add_axes([0.2, 0.2, 0.6, 0.6])
            swarmbox(x='day', y='score',
                     connect=True, connect_on=['pubid'],
                     data=plot_df.loc[ind],
                     order=['0', '3', '56', '59', '63', '112', '168'],
                     box_palette=['lightgray'],
                     swarm_color='k')
            plt.title('\n'.join(trt) + '\n', size=10)
            plt.xlabel('Study day')
            plt.ylabel(f'{modnames[m].upper()} Module Score')  
            pdf.savefig(figh)
            plt.close(figh)

with PngPdfPages(opj(out_folder, f'module_sex_pubid_lineplots.pdf')) as pdf:
    """Prep scores for summary plot with multiple modules. Mean across genes."""
    for trt in [trt34] + [[t] for t in treatments]:
        for m in modules:
            ind = plot_df['Treatment_Group'].isin(trt) & (plot_df['module'] == m)
            figh = plt.figure(figsize=(6, 4))
            swarmbox(x='day', y='score',
                     hue='sex',
                     hue_order=['female', 'male'],
                     data=plot_df.loc[ind],
                     order=['0', '3', '56', '59', '63', '112', '168'],
                     swarm_size=3)
            plt.title('\n'.join(trt))
            plt.xlabel('Participant sex')
            plt.ylabel(f'{modnames[m]} Module Score')  
            pdf.savefig(figh)
            plt.close(figh)

with PngPdfPages(opj(out_folder, f'module_sex2_pubid_lineplots.pdf')) as pdf:
    """Prep scores for summary plot with multiple modules. Mean across genes."""
    for trt in [trt34] + [[t] for t in treatments]:
        for m in modules:
            ind = plot_df['Treatment_Group'].isin(trt) & (plot_df['module'] == m)
            figh = plt.figure(figsize=(6, 4))
            swarmbox(x='sex', y='score',
                     hue='day',
                     hue_order=['0', '3', '56', '59', '63', '112', '168'],
                     connect=True, connect_on=['pubid'],
                     data=plot_df.loc[ind],
                     order=['female', 'male'],
                     swarm_size=3)
            plt.title('\n'.join(trt))
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Day')
            plt.xlabel('Participant sex')
            plt.ylabel(f'{modnames[m]} Module Score')  
            pdf.savefig(figh)
            plt.close(figh)


with PngPdfPages(opj(out_folder, f'module_lineplots.pdf')) as pdf:
    """Prep scores for summary plot with multiple modules.
    Mean across genes. Baseline subtract for each ptid. Mean across groups (3, 4)"""
    plot_df = pd.merge(samples, avg, how='left', on='sampleid')
    plot_df = plot_df.assign(Day=plot_df['day'].map(day_map))
    ind = plot_df['Treatment_Group'].isin(trt34)
    #ind = plot_df['Treatment_Group'] == 'Placebo'
    plot_df = plot_df.loc[ind].groupby(['module', 'pubid', 'Day'])['score'].agg('mean').reset_index()
    plot_df = pd.merge(plot_df, plot_df.loc[plot_df['Day'] == 0], how='left', on=['module', 'pubid'], suffixes=('', '_0'))
    plot_df = plot_df.assign(delta=plot_df['score'] - plot_df['score_0'])
    plot_df = plot_df.groupby(['module', 'Day'])['delta'].agg('mean').unstack('Day')
    
    """Multi-module plot, in different combinations of groups"""
    yl = (-0.5, 0.7)
    for mg in mod_groups:
        figh = plt.figure(figsize=(8, 4))
        axh = figh.add_axes([0.1, 0.1, 0.7, 0.8])
        for mod in mg:
            plt.plot(plot_df.columns, plot_df.loc[mod], '-s', label=mod, color=mod)
        xl = plt.xlim()
        plt.plot(xl, [0, 0], '--', color='gray')
        plt.xlim(xl)
        plt.ylim(yl)
        plt.xticks(ticks=day_ticks, labels=day_labels)
        plt.ylabel(f'Module (normalized counts)')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Module')
        pdf.savefig(figh)
        plt.close(figh)

    """Line plots of one module at a time, combining or separating groups"""
    for mod in modules:
        plot_df = pd.merge(samples, avg.loc[avg['module'] == mod], how='left', on='sampleid')
        plot_df = plot_df.assign(Day=plot_df['day'].map(day_map))
        
        figh = plt.figure(figsize=(8, 4))
        axh = figh.add_axes([0.1, 0.1, 0.5, 0.8])
        for trt in treatments:
            ind = plot_df['Treatment_Group'] == trt
            tmp = plot_df.loc[ind].set_index(['pubid', 'Day'])['score'].unstack('Day').mean(axis=0)
            plt.plot(tmp.index, tmp, '-s', label=trt)
        plt.xticks(ticks=day_ticks, labels=day_labels)
        plt.ylabel(f'{modnames[mod]} Module (normalized counts)')
        plt.xlabel('Study Day')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        pdf.savefig(figh)
        plt.close(figh)

        figh = plt.figure(figsize=(8, 4))
        axh = figh.add_axes([0.1, 0.1, 0.5, 0.8])
        for trt in trt34:
            ind = plot_df['Treatment_Group'] == trt
            tmp = plot_df.loc[ind].set_index(['pubid', 'Day'])['score'].unstack('Day').mean(axis=0)
            plt.plot(tmp.index, tmp, '-s', label=trt)
        plt.xticks(ticks=day_ticks, labels=day_labels)
        plt.xlabel('Study Day')
        plt.ylabel(f'{modnames[mod]} Module (normalized counts)')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        pdf.savefig(figh)
        plt.close(figh)

        figh = plt.figure(figsize=(8, 4))
        axh = figh.add_axes([0.1, 0.1, 0.5, 0.8])
        ind = plot_df['Treatment_Group'].isin(trt34)
        tmp = plot_df.loc[ind].set_index(['pubid', 'Day'])['score'].unstack('Day').mean(axis=0)
        plt.plot(tmp.index, tmp, '-s', color='black')
        plt.xticks(ticks=day_ticks, labels=day_labels)
        plt.xlabel('Study Day')
        plt.ylabel(f'{modnames[mod]} Module (normalized counts)')
        # plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        pdf.savefig(figh)
        plt.close(figh)





