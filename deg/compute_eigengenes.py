import pandas as pd
import numpy as np
import itertools
from functools import partial
import sys
import os
from scipy import stats

import sklearn

from os.path import join as opj

from fg_shared import *

"""SUMMARY:
Create eigengene scores with correct signs.

We decided this AM that we would use the eigengene approach. Summary argument is
that the PCA will weight genes that maximize the SHARED variance which will
tend to be the biological variation and not the technical variation.

It would be good to save that first PC1 vector and color the network plot nodes
for each module using these values weights for each gene. It should show that
the hub genes are the ones that are upweighted in the module score/PC1.

"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned_longform.csv')
out_folder = opj(project_folder, 'wgcna')
cts_fn = opj(project_folder, 'log_normalized_counts.csv')
# res_fn = opj(project_folder, 'agregated_results_2023-MAR-15.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

cts_df = pd.read_csv(cts_fn)
rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)

treatments = ['2 µg ID93 + 2 µg GLA-SE',
              '10 µg ID93 + 2 µg GLA-SE',
              '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
              '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)',
              'Placebo']

trt34 = ['2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)']

cts = cts_df.set_index('Unnamed: 0')
samples = pd.DataFrame(cts.columns, columns=['sampleid'])
samples = samples.assign(ptid=samples['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                         day=samples['sampleid'].map(lambda s: s.split('_')[-1]))
samples = pd.merge(samples, rx_df, how='left', on='ptid')

modules = modules_df['module'].unique()

def simple_avg(cts, genes):
    tmp = cts.loc[genes].mean(axis=0).T.reset_index()
    tmp.columns = ['sampleid', 'avg']
    return tmp

def eigengene_avg(cts, genes, training_ind):
    cts_a = cts.loc[genes, training_ind].T
    cts_b = cts.loc[genes, :].T

    """
    (1) Take subset of samples from treatments 3,4 and call these samples A. the whole 305 samples will be samples B
    (2) Compute mean and standard deviation for each gene among samples A. save these values as "A factors"
    (3) Subtract the mean and standard deviation for samples A using the A factors.
        Check that these now have precisely a mean of 0 and std of 1 for each gene
    (4) Compute PCA on just the genes from module 1 (repeat for all 8 modules).
        Save the first eigen vector (should be a one dimensional vector with one element per gene in the module)
    (5) Now go back to samples B and subtract the mean and divide by the std form "A factors".
        Check that mean for each gene is close to 0 but not exact and that std is close to 1 but not exactly
    (6) Then perform weighted average across genes in module 1 using the weights from the saved first eigen vector.
        You should now have one score per sample that is the eigengene for module 1
    (7) Save the "A factors" mean and std for each gene in a matrix and additionally save the
        weight that was applied to each gene from the first eigen vector.
        Save this matrix as a CSV to be able to plot weights and rescale gene appropriately later.
    """

    """a_means = cts_a.mean(axis=0)
    a_stds = cts_a.std(axis=0)

    cts_a_trans = (cts_a - a_means) / a_stds
    cts_b_trans = (cts_b - a_means) / a_stds"""

    scaler = sklearn.preprocessing.StandardScaler()
    scaler.fit(cts_a)
    cts_a_trans = scaler.transform(cts_a)

    pca = sklearn.decomposition.PCA(n_components=3, svd_solver='full', random_state=1)
    # pca.fit(cts_a_trans)
    check_a = pca.fit_transform(cts_a_trans)

    pca_b = pca.transform(scaler.transform(cts_b))
    pca_a = pca.transform(scaler.transform(cts_a))

    # np.isclose(check_a, pca_a)) TRUE

    eigengenes = pd.Series(pca_b[:, 0], index=cts.columns, name='eigengene')
    weights = pd.Series(pca.components_[0, :], index=genes, name='weight')

    return eigengenes, weights

avg = []
weights = []
for mod in modules:
    genes = modules_df.loc[modules_df['module'] == mod, 'gene'].tolist()
    tmp = simple_avg(cts, genes)
    tmp = tmp.assign(module=mod)

    training_ind = samples.loc[samples['Treatment_Group'].isin(trt34), 'sampleid']
    eg, w = eigengene_avg(cts, genes, training_ind)
    tmp = pd.merge(tmp, eg, left_on='sampleid', right_index=True)

    if tmp[['avg', 'eigengene']].corr().values[0, 1] < 0:
        """Flip sign when there is neg correlation of eigengene with avg"""
        tmp = tmp.assign(eigengene=-1*tmp['eigengene'])
        w = -1 * w
    tmp_scaled = sklearn.preprocessing.StandardScaler().fit_transform(tmp[['avg', 'eigengene']])
    tmp = tmp.assign(avg=tmp_scaled[:, 0],
                     eigengene=tmp_scaled[:, 1])
    avg.append(tmp)
    weights.append(w)
avg = pd.concat(avg, axis=0)
out_df = pd.merge(avg, samples, how='left', on='sampleid')

weights = pd.concat(weights, axis=0)

weights = pd.merge(modules_df[['module', 'gene']],
                   weights.reset_index().rename({'index':'gene'}, axis=1), on='gene')

# avg.groupby('module').apply(lambda gby:  gby[['avg', 'eigengene']].corr(method='spearman').values[0, 1])

# sns.jointplot(data=avg, x='avg', y='eigengene', hue='module')

# sns.histplot(data=weights, x='weight', hue='module')

weights.to_csv(opj(out_folder, 'modules_longform_weights.csv'))
out_df.to_csv(opj(out_folder, 'module_norm_cts.csv'))


