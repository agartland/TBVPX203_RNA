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
from myvolcano import plot_volcano

sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 12

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
res_fn = opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_res.csv')


res_df = pd.read_csv(res_fn)

sig_df = res_df.loc[(res_df['variable'] == 'day') & (res_df['cohort'] == '3,4') & res_df['sig']]
degs = sig_df['gene'].unique().tolist()
all_genes = res_df['gene'].unique().tolist()

"""Check first for genes associated with sex"""
plotdf = res_df.query("model=='~day + sex' & variable=='sex' & cohort=='3,4'") #" & (log2FC>0.5 | log2FC<-0.5) & pvalue<0.05 & FDR<0.2")
for comp, gby in plotdf.groupby('Comparison'):
    print(f"Significant genes for SEX at {comp}: {gby['sig'].sum()} of {gby['gene'].shape[0]}")

plotdf = res_df.query("model=='~day + sex' & variable=='day' & cohort=='3,4'")
plotdf = plotdf.assign(**{'Fold-change':plotdf['log2FC'].map(lambda v: 2**v)})

with PngPdfPages(opj(project_folder, f'volcanos.pdf')) as pdf:
    for comp, gby in plotdf.groupby('Comparison'):
        print(f"Significant genes for {comp}: {gby['sig'].sum()} of {gby['gene'].shape[0]} ({gby.loc[gby['Direction']=='UP', 'sig'].sum()} UP and {gby.loc[gby['Direction']=='DOWN', 'sig'].sum()} DOWN)")
        figh = plot_volcano(gby.rename({'FDR':'FDR-adj q-value'}, axis=1),
                             pvalue_col='FDR-adj q-value',
                             or_col='Fold-change',
                             sig_col='sig',
                             ann_col='gene',
                             annotate=None,
                             figsize=(5, 5),
                             fc_ticks=[1, 1.5, 2, 3, 4, 5])
        plt.xlim((1/6, 6))
        plt.ylim((1, 1e-20))
        plt.title(comp)
        pdf.savefig(figh)
        plt.close(figh)


"""
Significant genes for Day 112 vs. 56: 3 of (13687,) (3 UP and 0 DOWN)
Significant genes for Day 168 vs. 56: 0 of (13673,) (0 UP and 0 DOWN)
Significant genes for Day 3 vs. 0: 12 of (13902,) (7 UP and 5 DOWN)
Significant genes for Day 56 vs. 0: 2 of (13887,) (0 UP and 2 DOWN)
Significant genes for Day 59 vs. 56: 122 of (13817,) (89 UP and 33 DOWN)
Significant genes for Day 63 vs. 56: 313 of (13804,) (145 UP and 168 DOWN)
"""