import pandas as pd
import numpy as np

from fg_shared import _fg_data

mods = pd.read_csv(opj(_fg_data, "SCRI/TBVPX-203/RNA/21Nov2023/wgcna/new_wgcna_module_aligned_longform.csv"))
degs = pd.read_csv(opj(_fg_data, "SCRI/TBVPX-203/RNA/21Nov2023/degs_for_ingenuity/DEGs_day_sex_lme.csv"))

merged = pd.merge(degs, mods[['module', 'gene']], how='left', on='gene')

merged.to_csv(opj(_fg_data, "SCRI/TBVPX-203/RNA/21Nov2023/degs_for_ingenuity/DEGs_to_pathways_day_sex_lme.csv"))

print(merged.groupby(['module', 'Comparison2', 'Direction'])['gene'].count())

for i, gby in merged.groupby(['module', 'comparison']):
    fn = f'{i[0]}_{i[1]}_degs.csv'
    if gby.shape[0] >=9:
        print(fn, gby.shape[0])
        gby.to_csv(opj(_fg_data, f"SCRI/TBVPX-203/RNA/21Nov2023/degs_for_ingenuity/{fn}"))

    