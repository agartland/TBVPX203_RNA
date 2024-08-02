import pandas as pd
import numpy as np
import itertools
from functools import partial
import sys
import os
from scipy import stats
import matplotlib.pyplot as plt
from os.path import join as opj

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from adjustwithin import adjustwithin
from adjustwithin import adjustnonnan
from corrplots import partialcorr
from myboxplot import swarmbox


"""
Does MTB300 CD4+ T cell response at Day 0 correlate with gene modules at same time point?
"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'modules_longform_weights.csv')
out_folder = opj(project_folder, 'correlation_results')
scores_fn = opj(project_folder, 'wgcna', 'module_norm_cts.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

adata_folder = opj(_fg_data, 'SCRI/TBVPX-203/data/immune_adata')
ics_fn = opj(adata_folder, 'ics_analysis_23.csv')

rx_df = pd.read_csv(rx_fn)


modnames = {'turquoise':'COREG', 'brown':'BCELL', 'blue':'MITOSIS',
            'pink':'EOS', 'yellow':'IFN-I', 'green':'NEUTRO-I',
            'red':'NEUTRO-II', 'black':'MONO'}
modcolors = {v:k for k,v in modnames.items()}

ics = pd.read_csv(ics_fn, encoding='utf-8', encoding_errors='replace')

ind_cols = ['ptid', 'visitname', 'Treatment', 'tcellsub', 'antigen']
val_col = 'pctpos_adj'

ics = ics.assign(ptid=ics['ptid'].map(lambda v: f'{v:1.0f}'[:3] + '_' + f'{v:1.0f}'[3:]),
                 visitname=ics['visitno'].map({70:'Day 70', 56:'Day 56', 224:'Day 224', 0:'Day 0', 14:'Day 14', 28:'Day 28'}))

keep_days = ['Day 0', 'Day 14', 'Day 28', 'Day 56', 'Day 70', 'Day 224']
ics = ics.loc[ics['visitname'].isin(keep_days) & ics['antigen'].isin(['ID93', 'MTB 300'])]
ics = ics[ind_cols + [val_col]]


"""Load and prep module scores and compute difference variables"""
scores = pd.read_csv(scores_fn)

ind_cols = ['ptid', 'Treatment', 'module', 'day']
val_col = 'eigengene'

scores_0 = scores.loc[scores['day'] == 0]
scores_0 = scores_0.assign(sex_female=scores_0['sex'].map({'female':1, 'male':0}))

res = []
for module, gby in scores_0.groupby('module'):
    for ag in ['ID93', 'MTB 300']:
        ics_ss = ics.query(f'visitname=="Day 0" & tcellsub=="CD4+" & antigen=="{ag}"')
        tmp = pd.merge(gby[['ptid', 'eigengene', 'sex_female']],
                        ics_ss[['ptid', 'pctpos_adj']],
                        how='inner',
                        on='ptid')
        tmp = tmp.dropna()
        rho, pvalue = stats.spearmanr(tmp['pctpos_adj'], tmp['eigengene'])
        part_rho, part_pvalue = partialcorr(tmp['pctpos_adj'],
                                              tmp['eigengene'],
                                              adjust=[tmp['sex_female']],
                                              method='spearman')

        res.append(dict(
                        Antigen=ag,
                        Module=module,
                        Module_n=gby.dropna().shape[0],
                        n=tmp.shape[0],
                        rho=rho,
                        pvalue=pvalue,
                        part_rho=part_rho,
                        part_pvalue=part_pvalue))
res = pd.DataFrame(res)



"""Do MTB300 responses change over time?"""
ics_ss = ics.query('tcellsub=="CD4+" & antigen=="MTB 300"')
swarmbox(x='visitname', y='pctpos_adj', connect_on=['ptid'], connect=True, data=ics_ss.loc[ics_ss['visitname'].isin(['Day 0', 'Day 224'])])
tmp = ics_ss.set_index(['ptid', 'visitname', 'antigen'])['pctpos_adj'].unstack('visitname')
tmp = tmp[['Day 0', 'Day 224']].dropna()
wstat, pvalue = stats.wilcoxon(tmp['Day 0'], tmp['Day 224'])

"""Do module responses change comparing day 0 to day 168?"""
for module, gby in scores.groupby('module'):
    plt.figure()
    swarmbox(x='day', y='eigengene', connect_on=['ptid'], connect=True, data=gby.loc[gby['day'].isin([0, 168])])
    plt.title(module)
    tmp = gby.set_index(['ptid', 'day'])['eigengene'].unstack('day')
    tmp = tmp[[0, 168]].dropna()
    wstat, pvalue = stats.wilcoxon(tmp[0], tmp[168])
    print(f'{module}: p={pvalue:1.3g}')
