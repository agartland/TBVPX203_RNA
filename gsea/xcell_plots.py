import pandas as pd
import numpy as np
from os.path import join as opj
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

from fg_shared import *

sys.path.append(opj(_git, 'utils'))
from adjustwithin import adjustwithin, adjustnonnan
from myboxplot import swarmbox
from pngpdf import PngPdfPages

qdata_folder = opj(_trials, 'vaccine', 'IDRI', 'TBVPX203', 'qdata')
adata_folder = opj(_fg_data, 'SCRI', 'TBVPX-203', 'RNA')
results_folder = opj(adata_folder, '21Nov2023', 'deconvolution_results')

cts_fn = opj(adata_folder, '2019Dec', 'salmon_counts.csv')
rx_fn = opj(adata_folder, 'trt_pubid_2022-DEC-19.csv')
xcell_fn = opj(adata_folder, '21Nov2023', 'deconvolution_results', 'x_cell_results.csv')

rx = pd.read_csv(rx_fn)

cell_types = ['aDC',
             'B-cells',
             'Basophils',
             'CD4+ T-cells',
             'CD8+ T-cells',
             'cDC',
             'Class-switched memory B-cells',
             'DC',
             'Eosinophils',
             'iDC',
             'Memory B-cells',
             'MEP',
             'Monocytes',
             'naive B-cells',
             'Neutrophils',
             'NK cells',
             'NKT',
             'pDC',
             'Plasma cells',
             'pro B-cells',
             'Th1 cells']

plot_cells = ['Class-switched memory B-cells',
             'Eosinophils',
             'Monocytes',
             'Neutrophils',
             'Plasma cells',
             'Th1 cells']

xcell = pd.read_csv(xcell_fn)
xcell = xcell.set_index('Unnamed: 0').T.reset_index()
xcell = xcell.assign(Day=xcell['index'].map(lambda s: int(s.split('_')[-1])),
                     ptid=xcell['index'].map(lambda s: '_'.join(s.split('_')[1:3])))

xcell = pd.merge(xcell, rx[['ptid', 'pubid', 'Treatment_Group', 'Treatment']], how='left', on='ptid')

xcell = xcell.set_index(['pubid', 'Treatment_Group', 'Treatment', 'Day'])[cell_types].stack().reset_index().rename({'level_4':'Cell type', 0:'Proportion'}, axis=1)

comparisons = [(0, 3),
               (56, 59),
               (56, 63)]
test_df = xcell.loc[xcell['Treatment'].isin([3, 4])].set_index(['pubid', 'Treatment', 'Day', 'Cell type'])['Proportion'].unstack('Day').reset_index()

res = []
for cell in cell_types:
    for t1, t2 in comparisons:
        tmp = test_df.loc[test_df['Cell type'] == cell]
        n_left = tmp[t1].dropna().shape[0]
        n_right = tmp[t2].dropna().shape[0]
        
        tmp = tmp.dropna()
        n = tmp.shape[0]

        try:
            s, pvalue = stats.wilcoxon(tmp[t1].values, tmp[t2].values)
        except ValueError:
            s, pvalue = np.nan, np.nan
        res.append(dict(cell_type=cell,
                        left=t1,
                        right=t2,
                        n_left=n_left,
                        n_right=n_right,
                        n=n,
                        stat=s,
                        pvalue=pvalue,
                        avg_left=tmp[t1].mean(),
                        avg_right=tmp[t2].mean()))
res = pd.DataFrame(res)
res = res.assign(FDRq=adjustnonnan(res['pvalue'], method='fdr_bh'),
                 FWERp=adjustnonnan(res['pvalue'], method='holm'))

res.to_csv(opj(results_folder, 'xcell_fc_stats.csv'))


plot_df = test_df.loc[test_df['Cell type'].isin(plot_cells), ['pubid', 'Cell type', 0, 3, 56, 59, 63]]
# swarmbox(x='Cell type', hue='Day', y='Proportion', connect_on='pubid', connect=True, data=cib.loc[cib['Treatment'].isin([3, 4]) & cib['Cell type'].isin(['Neutrophils', 'Plasma cells'])])
with PngPdfPages(opj(results_folder, f'xcell_plots.pdf')) as pdf:
    for ct in plot_cells:
        figh = plt.figure(figsize=(5, 3))
        axh = figh.add_axes([0.2, 0.2, 0.6, 0.6])
        plt.plot(plot_df.loc[plot_df['Cell type'].isin([ct]), [0, 3, 56, 59, 63]].values.T, '-', lw=0.5, color='gray')
        plt.plot([0, 1, 2, 3, 4], plot_df.loc[plot_df['Cell type'].isin([ct]), [0, 3, 56, 59, 63]].mean(axis=0).values, '-s', lw=2, color='purple')
        plt.xticks([0, 1, 2, 3, 4], [0, 3, 56, 59, 63])
        plt.ylabel('Enrichment score')
        plt.xlabel('Study day')
        plt.title(ct, size='large')
        pdf.savefig(figh)
        plt.close(figh)

plot_df = (plot_df.set_index(['pubid', 'Cell type']) - plot_df[[56]].values).reset_index()
with PngPdfPages(opj(results_folder, f'xcell_fc_plots.pdf')) as pdf:
    for ct in plot_cells:
        figh = plt.figure(figsize=(5, 3))
        axh = figh.add_axes([0.2, 0.2, 0.6, 0.6])
        plt.plot(plot_df.loc[plot_df['Cell type'].isin([ct]), [56, 59, 63]].values.T, '-', lw=0.5, color='gray')
        plt.plot([0, 1, 2], plot_df.loc[plot_df['Cell type'].isin([ct]), [56, 59, 63]].mean(axis=0).values, '-s', lw=2, color='purple')
        plt.plot([0, 2], [0, 0], '--k', lw=1)
        plt.xticks([0, 1, 2], [56, 59, 63])
        plt.ylabel('Fold-change')
        plt.xlabel('Study day')
        plt.title(ct, size='large')
        pdf.savefig(figh)
        plt.close(figh)