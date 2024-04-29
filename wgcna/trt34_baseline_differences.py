import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools
from functools import partial
import sys
import os
from scipy import stats
import seaborn as sns
import statsmodels.formula.api as smf
from os.path import join as opj

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from myboxplot import swarmbox

"""No evidence that modules scores are significantly different at day 0 or day 56.
I also tried these models with sex interactively and adj for sex didn't change the result"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')

adata_folder = opj(_fg_data, 'SCRI/TBVPX-203/data/immune_adata')

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

delta = pd.read_csv(opj(adata_folder, 'module_adata.csv'))

scores_0 = delta.loc[delta['comparison'] == 'Day 0']
scores_56 = delta.loc[delta['comparison'] == 'Day 56']

swarmbox(x='module', y='delta', hue='Treatment', data=scores_56)

swarmbox(x='module', y='delta', hue='Treatment', data=scores_0)


tmp = scores_0.loc[scores_0['Treatment'].isin([3, 4])].set_index(['module', 'Treatment', 'ptid'])['delta'].unstack('module').reset_index()
tmp = tmp.assign(trt3=tmp['Treatment'] == 3)

for mod in modules_df.module.unique():
    print(smf.ols(formula=f'{mod} ~  trt3', data=tmp).fit().summary())


tmp = scores_56.loc[scores_56['Treatment'].isin([3, 4])].set_index(['module', 'Treatment', 'ptid'])['delta'].unstack('module').reset_index()
tmp = tmp.assign(trt3=tmp['Treatment'] == 3)

for mod in modules_df.module.unique():
    print(smf.ols(formula=f'{mod} ~ trt3', data=tmp).fit().summary())


"""ALSO QUICK PLOTS OF SEX DIFFERENCES FOR ICS< WBICS AND ELISA"""

tmp = pd.merge(ics.loc[(ics['Treatment'] == 3) & (ics['antigen'] == 'ID93') & (ics['tcellsub'] == 'CD4+')], rx_df[['sex', 'ptid']], how='left', on='ptid')
swarmbox(x='visitname', y='pctpos_adj', hue='sex', connect=True, connect_on=['ptid'], data=tmp)

tmp = pd.merge(wbics.loc[(wbics['Treatment'] == 3) & (wbics['antigen'] == 'ID93') & (wbics['tcellsub'] == 'CD4+')], rx_df[['sex', 'ptid']], how='left', on='ptid')
swarmbox(x='visitname', y='pctpos_adj', hue='sex', connect=True, connect_on=['ptid'], data=tmp)

tmp = pd.merge(elisa.loc[(elisa['rx_code'] == 'T3') & (elisa['analyte'] == 'Total IgG')], rx_df[['sex', 'ptid']], how='left', on='ptid')
swarmbox(x='visitname', y='MEPT', hue='sex', connect=True, connect_on=['ptid'], data=tmp)