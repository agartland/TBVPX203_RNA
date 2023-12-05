import pandas as pd
import numpy as np
import sys
import os

from os.path import join as opj

from fg_shared import *

"""
DESCRIPTION
Load raw-uncorrected salmon counts data: merge from two separate count files. Remove control samples.
Check for sample swaps in the HLA data.
Perform sample-swap, identify sample sex.
Save corrected gene count CSV
Save ptid-sex CSV
"""

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')

hla_fn = opj(project_folder, 'hla_genotypes.tsv')
in_cts_fn = opj(project_folder, 'salmon', 'salmon.merged.gene_counts.tsv')
in_tpm_fn = opj(project_folder, 'salmon', 'salmon.merged.gene_tpm.tsv')

"""Additional samples from a initial pilot batch that was successful but not repeated: should be merged"""
in_cts1_fn = opj(project_folder, 'salmon_23Apr2019', 'salmon.merged.gene_counts.tsv')
in_tpm1_fn = opj(project_folder, 'salmon_23Apr2019', 'salmon.merged.gene_tpm.tsv')

out_cts_fn = opj(project_folder, 'salmon', 'salmon.merged.gene_counts_corrected.tsv')
out_tpm_fn = opj(project_folder, 'salmon', 'salmon.merged.gene_tpm_corrected.tsv')
rx_fn = opj(project_folder, 'trt_pubid_2022-DEC-19.csv')
new_rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

def _hla2str(r):
    hstr = '_'.join(sorted([h[:7] for h in r[['A1', 'A2', 'B1', 'B2', 'C1', 'C2']]]))
    return hstr
hla = pd.read_csv(hla_fn, sep='\t')
hla = hla.loc[hla['subject'].str.slice(0, 1) != 'N']
hla = hla.assign(hla=hla.apply(_hla2str, axis=1),
                 ptid=hla['subject'].map(lambda s: '_'.join(s.split('_')[1:-1])))

"""print(hla.groupby(['ptid', 'hla'])['subject'].count())
Output shows that there is one clear sample swap"""

"""sampleid P203_903_3006_0 was swapped with P203_903_3005_0
Confirmed based on HLA genotyping data"""
swap_map = {'X203_903_3006_0':'X203_903_3005_0',
            'X203_903_3005_0':'X203_903_3006_0'}

rx = pd.read_csv(rx_fn)
cts = pd.read_csv(in_cts_fn, sep='\t')
tpm = pd.read_csv(in_tpm_fn, sep='\t')

sampleids = [c for c in cts.columns if c[0] == 'X']

"""Remove controls, swap samples"""
cts = cts.set_index('gene_id')[sampleids].rename(swap_map, axis=1)
tpm = tpm.set_index('gene_id')[sampleids].rename(swap_map, axis=1)

"""MERGE ADDITIONAL BATCH HERE"""
cts1 = pd.read_csv(in_cts1_fn, sep='\t')
tpm1 = pd.read_csv(in_tpm1_fn, sep='\t')
sampleids1 = [c for c in cts1.columns if c[0] == 'X']

cts1 = cts1.set_index('gene_id')[sampleids1]
tpm1 = tpm1.set_index('gene_id')[sampleids1]

cts = pd.merge(cts, cts1, how='left', on='gene_id').fillna(0)
tpm = pd.merge(tpm, tpm1, how='left', on='gene_id').fillna(0)

"""Output the de-swapped/corrected files"""
cts.reset_index().to_csv(out_cts_fn, index=False)
tpm.reset_index().to_csv(out_tpm_fn, index=False)


"""Identify participant sex"""
ygenes = ['DDX3Y', 'UTY', 'ZFY', 'PRKY']
"""These show a clear separation at logTPM = 0"""
# sns.pairplot(np.log(tpm.loc[['DDX3Y', 'UTY', 'ZFY', 'PRKY']].T + 0.1))
# sns.histplot(np.log(tpm.loc[['DDX3Y', 'UTY', 'ZFY', 'PRKY']].T.sum(axis=1)), bins=100)

sex = np.log(tpm.loc[ygenes].T.sum(axis=1) + 0.1).reset_index()
sex = sex.assign(sex=sex[0].map(lambda v: 'male' if v > 0 else 'female'))
sex = sex.rename({'index':'sampleid', 0:'log_tpm_y'}, axis=1)
sex = sex.assign(ptid=sex['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                 day=sex['sampleid'].map(lambda s: s.split('_')[-1]))
sex = sex[['ptid', 'sex']].drop_duplicates()

nrx = pd.merge(rx, sex, how='left', on='ptid')
nrx.to_csv(new_rx_fn, index=False)


