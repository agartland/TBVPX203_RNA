import pandas as pd
import sys
from os.path import join as opj
from fg_shared import _fg_data
from glob import glob

results_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023/gsea')

"""example file name: em_C7_red_vs_all.csv"""

res = []
for fn in glob(opj(results_folder, 'em_*.csv')):
    file = fn.split(os.path.sep)[-1]
    em, gset, module, vs, universe = file.split('_')
    tmp = pd.read_csv(fn)
    tmp = tmp.assign(gene_set=gset,
                     filename=file,
                     module=module,
                     universe=universe.replace('.csv', ''))
    if tmp.shape[0] > 0:
        res.append(tmp)
res = pd.concat(res, axis=0)

res.to_csv(opj(results_folder, 'aggregated_gsea_2024-MAR-14.csv'), index=False)

