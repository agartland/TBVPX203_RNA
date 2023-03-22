from fg_shared import *
from os.path import join as opj
import pandas as pd
import numpy as np
import sys
from scipy import stats

sys.path.append(opj(_git, 'utils'))
# from pngpdf import PngPdfPages
from embedding import *

sys.path.append(opj(_git, 'IDRI', 'TBVPX-203', 'RNA'))
from load_data import *

tpm = load_salmon(tpm=True)
rx = load_rx()

"""Focus on day 0 vs. day 3"""
tpm = tpm.loc[tpm['day'].isin([0, 3])]

"""Require that genes are counted >10 times in >50% of samples"""
summ = tpm.groupby('gene_id')['counts'].agg(lambda v: (v > 10).mean())
keepers = summ.index[(summ >= 0.5)]
tpm = tpm.loc[tpm['gene_id'].isin(keepers)]

"""Identify participants with both samples"""
tmp = tpm[['ptid', 'day']].drop_duplicates().groupby('ptid')['day'].count()
keep_ptids = tmp.index[tmp==2]
tpm = tpm.loc[tpm['ptid'].isin(keep_ptids)]

lfc = pd.merge(tpm.loc[tpm['day'] != 0], tpm.loc[tpm['day'] == 0][['counts', 'ptid', 'gene_id']], on=['ptid', 'gene_id'], suffixes=('', '_0'))
lfc = lfc.assign(lfc=np.log(lfc['counts'])-np.log(lfc['counts_0']))

plotdf = lfc.set_index(['gene_id', 'ptid', 'day'])['lfc'].unstack('gene_id').reset_index()
plotdf = pd.merge(plotdf, rx[['ptid', 'Treatment Group']], on='ptid', how='left').set_index(['ptid', 'day', 'Treatment Group'])

xy = embedObservations(plotdf, method='pca', n_components=2, metric='euclidean', downSample=None)
meta = xy.index.to_frame()

figh = plt.figure(1)
figh.clf()
axh = figh.add_subplot(plt.GridSpec(nrows=1, ncols=1, left=0.1, right=0.45)[0])

plotEmbedding(plotdf,
                  xyDf=xy,
                  labels=meta['day'],
                  plotLabels=False,
                  plotDims=[0, 1],
                  plotElipse=False,
                  method='PCA',
                  alpha=0.8,
                  plotLegend=True,
                  markerLabels=meta['Treatment Group'])