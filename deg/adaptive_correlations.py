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

"""TODO:
Make sure these correlations are done by first computing gene-wise differences (eg D56 from D59)
and then computing averages. Previously we used module eigen values which destroys this information.
"""

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from myboxplot import swarmbox

#sns.set_style('whitegrid')
mpl.rcParams['font.size'] = 12

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/2019Dec/Results')
modules_fn = opj(project_folder, 'updated_wgcna', 'WGCNA_modules_stricter_DEG_expression.csv')
out_folder = opj(project_folder, 'presentation_results')
cts_fn = opj(project_folder, 'updated_wgcna', 'normalized_data_all.csv')
res_fn = opj(project_folder, 'agregated_results_2023-MAR-15.csv')
rx_fn = opj(_fg_data, 'SCRI/TBVPX-203/RNA/trt_pubid_2022-DEC-19.csv')

res_df = pd.read_csv(res_fn)
cts_df = pd.read_csv(cts_fn)
rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)
modules = ['turquoise', 'brown', 'green', 'yellow', 'blue','red','black']
modules_df = modules_df[modules].stack().reset_index().rename({'level_1':'module', 0:'gene'}, axis=1)[['module', 'gene']]

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

comparison_map = {'0_3':'Day 3 vs. 0',
                  '56_59':'Day 59 vs. 56',
                  '56_63':'Day 63 vs. 56'}

cts = cts_df.set_index('Unnamed: 0')
samples = pd.DataFrame(cts.columns, columns=['sampleid'])
samples = samples.assign(ptid=samples['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                         day=samples['sampleid'].map(lambda s: s.split('_')[-1]))
samples = pd.merge(samples, rx_df, how='left', on='ptid')


day_ticks = [0, 1, 2, 3, 4, 5, 6]
day_labels = ['0', '3', '56', '59', '63', '112', '168']
day_map = {k:v for k,v in zip(day_labels, day_ticks)}

mod_groups = [('blue', 'black', 'brown', 'red','green', 'turquoise', 'yellow'),
              ('blue',  'brown', 'red'),
              ('green', ),
              ('yellow','black','turquoise')]

plot_df.to_csv(opj(project_folder, 'module_norm_cts_avg.csv'))