import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# import seaborn as sns
import itertools
import scipy.cluster.hierarchy as sch
import sys
import os

from os.path import join as opj
from glob import glob

import networkx as nx
import holoviews as hv
# import holoviews.operation.datashader
from holoviews import opts
# hv.extension('bokeh')
# from fa2 import ForceAtlas2
# from scipy import sparse

from fg_shared import *

"""SUMMARY:
Create network plots of DEGs, potentially separated by up/down regulation or WGCNA modules,
depending on complexity.
"""

sys.path.append(opj(_git, 'utils'))
from pngpdf import PngPdfPages
from cluster_alignment import align_clusters

#sns.set_style('whitegrid')

project_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/21Nov2023')
modules_fn = opj(project_folder, 'wgcna', 'new_wgcna_module_aligned_longform.csv')
out_folder = opj(project_folder, 'module_results')
# out_file = opj(out_folder, f'deg_networks.pdf')
out_file = opj(out_folder, f'deg_networks_pos_05.pdf')
cts_fn = opj(project_folder, 'log_normalized_counts.csv')
res_fn = opj(project_folder, 'degs_day_sex', 'agg_day_sex_lme_res.csv')
rx_fn = opj(project_folder, 'trt_pubid_2023-NOV-28.csv')

# res_df = pd.read_csv(res_fn)
cts_df = pd.read_csv(cts_fn)
rx_df = pd.read_csv(rx_fn)

modules_df = pd.read_csv(modules_fn)
"""Somehow there are three genes that made the filter for DEG testing
but missed the filter for normalized counts matrix: dropping them here but should be added back to normalized cts"""
modules_df = modules_df.loc[~modules_df['gene'].isin(['MYO5C','DDX11L5', 'LOC101928891'])]

# from matplotlib_venn import venn2
# venn2([set(omodules_df['gene']), set(modules_df['gene'])], set_labels=('Old modules', 'New modules'), set_colors=('r', 'b'))

treatments = ['2 µg ID93 + 2 µg GLA-SE', 'Placebo',
               '2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)',
               '10 µg ID93 + 2 µg GLA-SE']
trt34 = ['2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)',
               '2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)']

corr_threshold = 0.5

res_df = pd.read_csv(res_fn)
sig_df = res_df.loc[(res_df['variable'] == 'day') & (res_df['cohort'] == '3,4') & res_df['sig']]
sig_df = pd.merge(sig_df, modules_df, how='left', on='gene')
sig_df = sig_df.dropna(subset=['module'])
sig_genes = sig_df['gene'].unique().tolist()

cts = cts_df.set_index('Unnamed: 0').loc[sig_genes]

samples = pd.DataFrame(cts.columns, columns=['sampleid'])
samples = samples.assign(ptid=samples['sampleid'].map(lambda s: '_'.join(s.split('_')[1:-1])),
                         day=samples['sampleid'].map(lambda s: s.split('_')[-1]))
samples = pd.merge(samples, rx_df, how='left', on='ptid')

samps34 = samples.loc[samples['Treatment_Group'].isin(trt34)]

"""Mat has genes as features (columns) and observations as rows"""
mat = cts[samps34['sampleid']].T

"""Standardize so each feature (column) has mean = 0 and std = 1"""
smat = (mat - np.mean(mat, axis=0)) / np.std(mat, axis=0)

"""Empirical covariance matrix (426 genes x genes)"""
cov_smat = np.cov(smat.T)

"""Using method='pearson' here makes it identical to the empirical covariance matrix"""
corr_mat = mat.corr(method='spearman')

adjacency = corr_mat > corr_threshold
# adjacency = np.abs(corr_mat) > corr_threshold
adjacency.values[np.diag_indices_from(adjacency)] = 0

'''
from gglasso.problem import glasso_problem
from gglasso.helper.utils import sparsity, zero_replacement, normalize, log_transform
from gglasso.helper.basic_linalg import scale_array_by_diagonal

N = mat.shape[0]

P = glasso_problem(cov_smat, N, reg_params=dict(lambda1=0.05), latent=False, do_scaling=False)
print(P)

lambda1_range = np.logspace(0, -3, 30)
modelselect_params = dict(lambda1_range=lambda1_range)

P.model_selection(modelselect_params=modelselect_params,
                  method='eBIC',
                  gamma = 0.1)

# regularization parameters are set to the best ones found during model selection
print(P.reg_params)

# tmp = P.modelselect_stats
sol = P.solution.precision_
"""The precision matrix is the inverse of the covariance matrix,
where 0 means two variables are conditionally independent, after adjusting
for all other variables"""


S0 = np.cov(smat.T, bias=True)
"""Not needed becuase I pre-standardized the data"""
S = scale_array_by_diagonal(S0)

P = glasso_problem(S, N, latent=True, do_scaling=False)
print(P)

lambda1_range = np.logspace(0.5, -1.5, 8)
mu1_range = np.logspace(1.5, -0.2, 6)

modelselect_params = {'lambda1_range': lambda1_range, 'mu1_range': mu1_range}

P.model_selection(modelselect_params=modelselect_params, method='eBIC', gamma=0.25)

print(P.reg_params)

reg_params = dict(lambda1=0.22758459260747887, mu1=31.622776601683796)
reg_params = dict(lambda1=0.5, mu1=31.62)
P = glasso_problem(S, N, reg_params=reg_params, latent=True, do_scaling=False)
P.solve()
P.solution.calc_adjacency(t=1e-2)
adjacency = P.solution.adjacency_

print(adjacency.mean())

prec_est = P.solution.precision_
prec_est[np.diag_indices_from(prec_est)] = 0

'''

nodes_df = sig_df.groupby('gene').apply(lambda df: df.loc[df['log2FC'].idxmax()])
nodes_df = nodes_df.loc[corr_mat.index, ['Comparison', 'gene', 'module', 'log2FC', 'pvalue', 'FDR']]
nodes_df = nodes_df.assign(logP=nodes_df['pvalue'].map(np.log10),
                            Direction=nodes_df['log2FC'].map(lambda v: {True:'UP', False:'DOWN'}[v>0]),
                            Module=nodes_df['module'])
nodes_df = nodes_df.reset_index(drop=True)


direction_colors = {'UP':'red', 'DOWN':'blue'}
module_colors = {k:k for k in ['turquoise', 'brown', 'green', 'yellow', 'blue','red','black', 'grey', 'pink']}
comparison_colors = {'Day 3 vs. 0':'magenta',
                     'Day 59 vs. 56':'red',
                     'Day 63 vs. 56':'black',
                     'Day 112 vs. 56':'orange',
                     'Day 56 vs. 0':'brown'}

def _create_graph(adjacency, nodes_df, edge_weights, layout='neato'):
    g = nx.Graph(adjacency)

    """Setting as node attributes here is important for using in bokeh plots later"""
    nx.set_node_attributes(g, nodes_df.T.to_dict())
    # edge_attr = {(i, j):{'dist':edge_weights[i, j]} for i, j in g.edges}
    # nx.set_edge_attributes(g, edge_attr)

    pos = nx.nx_agraph.graphviz_layout(g, prog=layout)
    # pos = nx.fruchterman_reingold_layout(G_pub)
    # pos = forceatlas2.forceatlas2_networkx_layout(G_sub, iterations=200)
    return g, pos

def _render_hv(g, pos, color_col, color_map):
    hv_graph = hv.Graph.from_networkx(g, pos)
    """Render the graph using Holoviews and save the output to an HTML file"""
    hv_graph = hv_graph.opts(width=1200, height=1200, node_size=10, edge_line_width=0.5, edge_color='gray',
                             node_line_color='gray', node_color=color_col, cmap=color_map, colorbar=True)
    # pdf.savefig(hv.render(hv_graph, backend='matplotlib'))
    return hv_graph

def _render_mpl(g, pos, nodes_df, color_col, color_map=None, node_size=30, figsize=(10, 8), **draw_params):
    if color_map is None:
        color_map = {k:k for k in nodes_df[color_col].unique()}

    figh = plt.figure(figsize=figsize)
    axh = figh.add_axes([0.0, 0.0, 0.7, 0.9])
    if type(color_map) is dict:
        nx.draw(g, pos,
                node_color=nodes_df[color_col].map(color_map),
                node_size=[node_size for x in g.nodes],
                edge_color='gray',
                with_labels=False, **draw_params)

        key_cts = nodes_df.groupby(color_col).count().iloc[:, 0].sort_values(ascending=False)
        plt.legend([plt.Circle(xy=(0,0), color=color_map[k]) for k in key_cts.index],
                   [f'{k} (n={key_cts[k]})' for k in key_cts.index],
                   bbox_to_anchor=(1, 1), loc='upper left', title=color_col)
    else:
        nx.draw(g, pos,
                node_color=nodes_df[color_col],
                cmap=color_map,
                node_size=[node_size for x in g.nodes],
                edge_color='lightgray',
                with_labels=False, **draw_params)
    return figh

# renderer = hv.renderer('bokeh')
with PngPdfPages(out_file) as pdf:
    g, pos = _create_graph(adjacency, nodes_df, edge_weights=corr_mat.values)

    figh = _render_mpl(g, pos, nodes_df, color_col='Module', color_map=module_colors)
    pdf.savefig(figh)
    plt.close(figh)

    figh = _render_mpl(g, pos, nodes_df, color_col='Direction', color_map=direction_colors)
    pdf.savefig(figh)
    plt.close(figh)

    figh = _render_mpl(g, pos, nodes_df, color_col='Comparison', color_map=comparison_colors)
    pdf.savefig(figh)
    plt.close(figh)

    # hv_graph = _render_hv(g, pos, color_col='Module', color_map=module_colors)
    # renderer.save(hv_graph, opj(out_folder, f'deg_networks'))


    for mod in nodes_df['module'].unique():
        mod_genes = nodes_df['gene'].loc[nodes_df['module'] == mod].tolist()
        tmp_mat = corr_mat.loc[mod_genes, :].loc[:, mod_genes]
        """Prune the genes that have too few connections"""
        keep_genes = tmp_mat.loc[(tmp_mat > corr_threshold).sum(axis=1) >= 3].index
        # keep_genes = tmp_mat.loc[(np.abs(tmp_mat) > corr_threshold).sum(axis=1) >= 3].index
        tmp_mat = corr_mat.loc[keep_genes, :].loc[:, keep_genes]

        # adjacency = np.abs(tmp_mat.values) > corr_threshold
        adjacency = tmp_mat.values > corr_threshold
        """Adding weights didn't change the neato layout"""
        #adjacency = tmp_mat.values
        #adjacency[tmp_mat.values <= corr_threshold] = 0

        adjacency[np.diag_indices_from(adjacency)] = 0
        tmp_nodes = nodes_df.set_index('gene', drop=False).loc[tmp_mat.index]
        tmp_nodes = tmp_nodes.reset_index(drop=True)

        g, pos = _create_graph(adjacency, tmp_nodes, edge_weights=tmp_mat.values, layout='neato')
        mx = np.quantile(tmp_nodes['log2FC'].abs(), [0.9])[0]

        centrality = nx.betweenness_centrality(g, endpoints=True)
        top10 = tmp_nodes.loc[pd.Series(centrality).sort_values(ascending=False).index[:15]]
        label_inds = top10.index

        figh = _render_mpl(g, pos, tmp_nodes, color_col='log2FC', color_map=mpl.cm.PuOr_r, vmin=-mx, vmax=mx)
        plt.title(mod.title(), size=20)
        labs = []
        for i in label_inds:
            plt.annotate(tmp_nodes.loc[i, 'gene'], xy=pos[i], size=9)
            labs.append(tmp_nodes.loc[i, 'gene'])
        plt.annotate('\n'.join(labs), xy=(0, 1), va='top', ha='left', xycoords='axes fraction', size=12)    

        axh = figh.add_axes([0.7, 0.2, 0.05, 0.6])
        figh.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-mx, vmax=mx), cmap=mpl.cm.PuOr_r),
                      cax=axh, orientation='vertical', label='log2-FC')
        pdf.savefig(figh)
        plt.close(figh)


