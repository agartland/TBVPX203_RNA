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


"""
Longitudinal change in modules scores
DIDNT COMPLETE THIS. HS computed using lme4 models in R

"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'modules_longform_weights.csv')
out_folder = opj(project_folder, 'correlation_results')
scores_fn = opj(project_folder, 'wgcna', 'module_norm_cts.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

scores = pd.read_csv(scores_fn)

ind_cols = ['ptid', 'Treatment', 'module', 'day']
val_col = 'eigengene'

scores_0 = scores.loc[scores['day'] == 0]
scores_n0 = scores.loc[scores['day'] != 0]

delta_0 = pd.merge(scores_n0[ind_cols + [val_col]],
                scores_0[ind_cols[:3] + [val_col]],
                how='inner',
                on=ind_cols[:3],
                suffixes=('', '_0'))

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
