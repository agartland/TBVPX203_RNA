import pandas as pd
import numpy as np
import sys
from os.path import join as opj

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

sns.set_style('whitegrid')

from fg_shared import _fg_data

sys.path.append(opj(_git, 'utils'))
from myboxplot import swarmbox
from myvolcano import altair_scatter

results_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/2019Dec/Results')
res = pd.read_csv(opj(results_folder, 'agregated_results_2023-MAR-15.csv'))

def _compare_methods(res, methods, comparison, col, figh=None, pvalue=0.05, FDR=0.2, log2FC=0.5):
    ind = res['method'].isin(methods) & (res['comparison'] == comparison)
    tmp = res.loc[ind].set_index(['gene', 'method'])

    tmp = tmp.assign(sig=((tmp['pvalue'] < pvalue) & (tmp['FDR'] < FDR) & (np.abs(tmp['log2FC']) > log2FC)).astype(int) * np.sign(tmp['log2FC']))

    if col in ['pvalue', 'FDR']:
        NA = 1
        transform = 'log'
    elif col in ['tstat', 'beta', 'log2FC']:
        NA = 0
        transform = 'linear'

    plotdf = tmp[col].unstack('method').fillna(NA).reset_index()

    sig = tmp['sig'].unstack('method').fillna(0).reset_index()
    
    mapping = {(1, 1): 'Both up',
               (-1, 1): 'X down, Y up',
               (1, -1): 'X up, Y down',
               (-1, -1): 'Both down',
               (1, 0):'X-only up',
               (0, 1):'Y-only up',
               (-1, 0):'X-only down',
               (0, -1):'Y-only down',
               (0, 0):'Neither'}

    plotdf = plotdf.assign(sig=sig.apply(lambda r: mapping[(r[methods[0]], r[methods[1]])], axis=1))

    if figh is None:
        plotdf = plotdf.loc[plotdf['sig']!='Neither']
        figh = altair_scatter(methods[0], methods[1], 'sig', data=plotdf, tooltip=['gene'], yscale=transform, xscale=transform, palette=None, size=20, stroke=2, fontsize=14, title='', reversex=False, reversey=False)
    else:
        axh = figh.add_axes([0.15, 0.15, 0.75, 0.75], yscale=transform, xscale=transform)
        plt.grid(True, linewidth=1)
        axh.yaxis.set_minor_locator(AutoMinorLocator())

        for color, lab in zip(mpl.cm.tab10.colors, plotdf['sig'].unique()):
            if lab == 'Neither':
                z = -1
                alpha = 0.3
                s = 5
            else:
                z = 2
                alpha = 0.6
                s = 20
            plt.scatter(plotdf.loc[plotdf['sig'] == lab, methods[0]],
                        plotdf.loc[plotdf['sig'] == lab, methods[1]],
                        color=color, alpha=alpha, s=s, zorder=z,
                        label=lab)

        plt.ylabel(methods[1])
        plt.xlabel(methods[0])
        plt.legend()
    return figh, plotdf

figh, plotdf = _compare_methods(res, methods=['voom', 'voom_dup_corr'],
                         comparison='56_63',
                         col='log2FC',
                         figh=plt.figure(figsize=(8, 6)), pvalue=0.05, FDR=0.2, log2FC=0.5)
#plt.xlim((-2.5, 2.5))
#plt.ylim((-2.5, 2.5))

figh, plotdf = _compare_methods(res, methods=['voom', 'kimma'],
                         comparison='56_59',
                         col='pvalue',
                         figh=plt.figure(figsize=(8, 6)), pvalue=0.05, FDR=0.2, log2FC=0.5)

mpl.rcParams['font.size'] = 14

with PngPdfPages(opj(results_folder, 'compare_deg_methods.pdf')) as pdf:
    for comp in ['0_3', '56_59', '56_63']:
        for m1, m2 in itertools.combinations(['deseq', 'dream', 'deseq_day', 'kimma', 'voom', 'voom_dup_corr'], 2):
            figh, plotdf = _compare_methods(res, methods=[m1, m2],
                             comparison=comp,
                             col='log2FC',
                             figh=plt.figure(figsize=(8, 6)), pvalue=0.05, FDR=0.2, log2FC=0.5)
            plt.title(f'{m1} vs. {m2} for {comp} comparison')
            pdf.savefig(figh)
            figh, plotdf = _compare_methods(res, methods=[m1, m2],
                             comparison=comp,
                             col='pvalue',
                             figh=plt.figure(figsize=(8, 6)), pvalue=0.05, FDR=0.2, log2FC=0.5)
            plt.title(f'{m1} vs. {m2} for {comp} comparison')
            pdf.savefig(figh)