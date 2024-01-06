import pandas as pd
import numpy as np
import itertools
import sys
import os

from os.path import join as opj

from fg_shared import *

"""SUMMARY:
Re-assign colors names to the new WGCNA modules so they closely resemble the old modules.
"""

sys.path.append(opj(_git, 'utils'))
from cluster_alignment import align_clusters

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_THR6.csv')
realigned_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned.csv')
realigned_longform_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned_longform.csv')

"""Load the old module definitions"""
omodules_df = pd.read_csv(opj(_fg_data, 'SCRI/TBVPX-203/RNA/2019Dec/Results/WGCNA_modules_stricter_DEG_expression.csv'))
omodules = omodules_df.columns.tolist()
omodules_df = omodules_df[omodules].stack().reset_index().rename({'level_1':'module', 0:'gene'}, axis=1)[['module', 'gene']]

"""Load the new modules"""
modules_df = pd.read_csv(modules_fn)
modules = modules_df.columns.tolist()[1:] #exclude "Unnamed: 0" column
modules_df = modules_df[modules].stack().reset_index().rename({'level_1':'module', 0:'gene'}, axis=1)[['module', 'gene']]

"""Align on shared genes"""
shgenes = list(set(omodules_df['gene']).intersection(set(modules_df['gene'])))
sh_modules = modules_df.set_index('gene').loc[shgenes]
sh_omodules = omodules_df.set_index('gene').loc[shgenes]
remapped_sh_module_labels = align_clusters(sh_omodules['module'].values, sh_modules['module'].values)

new_to_old_mapper = {new:rm for rm, new in zip(remapped_sh_module_labels, sh_modules['module'].values)}

modules_df = modules_df.assign(module=modules_df['module'].map(new_to_old_mapper))
modules_df.to_csv(realigned_longform_fn, index=False)

realigned = pd.read_csv(modules_fn).drop(['Unnamed: 0'], axis=1).rename(lambda c: new_to_old_mapper.get(c, c), axis=1)
realigned.to_csv(realigned_fn, index=True)
