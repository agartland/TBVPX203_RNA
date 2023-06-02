import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import numpy as np
import itertools

#import holoviews as hv

from os.path import join as opj

from fg_shared import _fg_data, _git

genesets_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/nanostring_genesets')

from glob import glob

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages

files = {'TBVPX_degs_6K':'WGCNA_modules_DEG_expression.csv',
         'TBVPX_degs_369':'WGCNA_modules_stricter_DEG_expression.csv',
         'NSTB107':'NSTB107_Johnson.csv'}

genes = {}


probesets = []
for fn in glob(opj(genesets_folder, '*.xlsx')):
    name = fn.split(os.path.sep)[-1].replace('.xlsx', '')
    probesets.append(name)
    files.update({name:fn.split(os.path.sep)[-1]})

    tmp = pd.read_excel(fn)
    genes[name] = set(tmp['Gene'].values)


degs6k = pd.read_csv(opj(genesets_folder, files['TBVPX_degs_6K'])).iloc[:, 1:].stack().reset_index().rename({'level_1':'cluster', 0:'gene'}, axis=1).drop('level_0', axis=1)
degs = pd.read_csv(opj(genesets_folder, files['TBVPX_degs_369'])).iloc[:, 1:].stack().reset_index().rename({'level_1':'cluster', 0:'gene'}, axis=1).drop('level_0', axis=1)

genes['TBVPX_degs_6K'] = set(degs6k['gene'])
genes['TBVPX_degs_369'] = set(degs['gene'])

genes['NSTB107'] = set(pd.read_csv(opj(genesets_folder, files['NSTB107']), header=None).iloc[0, :].values)

ns_keys = ['NS_Hs_HostResponse Genes',
             'NS_Hs_Myeloid_V2.0 Genes',
             'Type I interferon signaling 713 Genes',
             'Immune response to viral infection 339 Genes',
             'Th17 response 696 Genes',
             'Lymphocyte activation 678 Genes',
             'NS_Immunology_v2_C2328 Genes',
             'Interferon signaling 525 Genes']

with PngPdfPages(opj(genesets_folder, 'venn_diagrams.pdf')) as pdf:
    for l in ns_keys:
        for tbvpx in ['TBVPX_degs_369', 'TBVPX_degs_6K']:
            labels = [tbvpx, 'NSTB107', l]
            figh = plt.figure()
            venn3(subsets=tuple(genes[g] for g in labels),
                  set_labels=labels, ax=figh.add_axes([0.1, 0.1, 0.8, 0.8]))
            pdf.savefig(figh)

"""

tmp = []
keys = list(genes.keys())
for i,j in itertools.combinations(range(len(keys)), 2):
    if len(genes[keys[j]]) > len(genes[keys[i]]):
        ii = j
        jj = i
    else:
        ii = i
        jj = j
    tmp.append(dict(source=keys[ii],
                    target=keys[jj],
                    value=len(genes[keys[ii]].intersection(genes[keys[jj]]))))
links = pd.DataFrame(tmp)


hv.extension('bokeh')
hv.output(size=200)

nodes = pd.DataFrame([{'index':i, 'name':keys[i]} for i in range(len(keys))])

nodes = hv.Dataset(nodes, 'index')

chord = hv.Chord((links, nodes)).select(value=(5, None))
chord.opts(
    opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(), 
               labels='name', node_color=dim('index').str()))
"""