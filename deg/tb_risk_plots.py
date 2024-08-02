import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools
import sys
import os
from scipy import stats
import seaborn as sns

from os.path import join as opj

from fg_shared import *

"""SUMMARY:
TB risk score boxplots
"""

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from myboxplot import swarmbox

#sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 12

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
scores_fn = opj(project_folder, 'tb_risk_signature', 'Tb_risk_signature_results.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')
out_folder = opj(project_folder, 'tb_risk_signature')

scores_df = pd.read_csv(scores_fn).rename({'Unnamed: 0':'sampleid'}, axis=1)
rx_df = pd.read_csv(rx_fn)

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

trt34 = ['2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)']

scores_df = scores_df.assign(ptid=scores_df['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                             day=scores_df['sampleid'].map(lambda s: s.split('_')[-1]))
scores_df = pd.merge(scores_df, rx_df, how='left', on='ptid')

day_ticks = [0, 1, 2, 3, 4, 5, 6]
day_labels = ['0', '3', '56', '59', '63', '112', '168']
day_map = {k:v for k,v in zip(day_labels, day_ticks)}

with PngPdfPages(opj(out_folder, f'tb_risk_boxplots.pdf')) as pdf:
    """Prep scores for summary plot with multiple modules. Mean across genes."""
    for trt in [trt34] + [[t] for t in treatments]:
        for m in ['Darboe_RISK_11', 'Sweeney_OD_3', 'Thompson_FAIL_13']:
            ind = scores_df['Treatment_Group'].isin(trt)
            figh = plt.figure(figsize=(6, 4))
            axh = figh.add_axes([0.2, 0.2, 0.6, 0.6])
            swarmbox(x='day', y=m,
                     connect=True, connect_on=['pubid'],
                     data=scores_df.loc[ind],
                     order=['0', '3', '56', '59', '63', '112', '168'],
                     box_palette=['lightgray'],
                     swarm_color='k',
                     axh=axh)
            plt.title('\n'.join(trt))
            plt.xlabel('Study day')
            plt.ylabel(f'{m} Score')
            pdf.savefig(figh)
            plt.close(figh)