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
Compare estimated log-FC for day 3 vs 0 and day 59 vs 56
"""

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages

#sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 12

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
out_folder = opj(project_folder, 'degs_day_sex')
res_fn = opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_res.csv')
sig_res_fn = opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_sig_ms.csv')

sig = pd.read_csv(sig_res_fn)
sig = sig.loc[(sig['variable'] == 'day') & sig['Comparison'].isin(['Day 3 vs. 0', 'Day 59 vs. 56'])]
degs = sig['gene'].unique().tolist()
degs_3 = sig.loc[sig['Comparison'].isin(['Day 3 vs. 0']), 'gene'].unique().tolist()
degs_59 = sig.loc[sig['Comparison'].isin(['Day 59 vs. 56']), 'gene'].unique().tolist()
degs_3only = list(set(degs_3).difference(set(degs_59)))
degs_59only = list(set(degs_59).difference(set(degs_3)))
degs_both = list(set(degs_59).intersection(set(degs_3)))

res = pd.read_csv(res_fn)
res = res.loc[(res['variable'] == 'day') & (res['model'] == '~day + sex') & res['gene'].isin(degs) & (res['cohort'] == '3,4')]

cols = ['gene', 'statistic', 'log2FC', 'pvalue', 'FDR', 'Comparison', 'logP']
jres = pd.merge(res.loc[(res['Comparison'] == 'Day 3 vs. 0'), cols],
                res.loc[(res['Comparison'] == 'Day 59 vs. 56'), cols],
                how='outer',
                on='gene',
                suffixes=('_3', '_59'))

with PngPdfPages(opj(out_folder, f'deg_day3_day_59.pdf'), dpi=200) as pdf:
    figh = plt.figure(figsize=(7, 4))
    params = dict(x='log2FC_3', y='log2FC_59', edgecolor='gray', s=25)
    plt.scatter(data=jres.loc[jres['gene'].isin(degs_59only)], color='lightgray', label='Day 59 only', **params)
    plt.scatter(data=jres.loc[jres['gene'].isin(degs_3only)], color='dodgerblue', label='Day 3 only', **params)
    plt.scatter(data=jres.loc[jres['gene'].isin(degs_both)], color='salmon', label='Day 3 and 59', **params)
    plt.plot([-2, 2], [-2, 2], '--', color='gray')
    plt.xlabel('Day 3 vs. 0\n(log2-FC)')
    plt.ylabel('Day 59 vs. 56\n(log2-FC)')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='DEG at')
    plt.tight_layout()
    plt.xlim((-2, 2))
    plt.ylim((-2, 2))
    pdf.savefig(figh)

    figh = plt.figure(figsize=(7, 4))
    params = dict(x='pvalue_3', y='pvalue_59', edgecolor='gray', s=25)
    plt.scatter(data=jres.loc[jres['gene'].isin(degs_59only)], color='lightgray', label='Day 59 only', **params)
    plt.scatter(data=jres.loc[jres['gene'].isin(degs_3only)], color='dodgerblue', label='Day 3 only', **params)
    plt.scatter(data=jres.loc[jres['gene'].isin(degs_both)], color='salmon', label='Day 3 and 59', **params)
    plt.plot((1e-17, 1), (1e-17, 1), '--', color='gray')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Day 3 vs. 0\n(p-value)')
    plt.ylabel('Day 59 vs. 56\n(p-value)')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='DEG at')
    plt.tight_layout()
    plt.xlim((1e-17, 1))
    plt.ylim((1e-17, 1))
    pdf.savefig(figh)