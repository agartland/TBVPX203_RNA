import pandas as pd
import sys
from os.path import join as opj
from fg_shared import _fg_data

results_folder = opj(_fg_data, 'SCRI/TBVPX-203/RNA/2019Dec/Results')

methods = dict(deseq='Deg_resultsdeseq',
               dream='Deg_resultsDream',
               deseq_day='deseq_day',
               kimma='kimma',
               voom='limma_voom',
               voom_dup_corr='limma_voom_duplicate_correlation')
day_comparisons = ['0_3', '56_59', '56_63']

files = dict(deseq={c:f'res_{c}_deseq.csv' for c in day_comparisons},
               dream={c:f'res_{c}_Dream.csv' for c in day_comparisons},
               deseq_day={c:f'res_{c}_deseq.csv' for c in day_comparisons},
               kimma={c:f'lme.contrast_{c}_kimma.csv' for c in day_comparisons},
               voom={c:f'res_{c}_limma_voom.csv' for c in day_comparisons},
               voom_dup_corr={c:f'res_{c}_limma_voom_duplicate_correlation.csv' for c in day_comparisons})

col_map = {'pval':'pvalue',
            'DEG_results.logFC':'log2FC',
            'DEG_results.AveExpr':'avg_expression',
            'DEG_results.t':'tstat',
            'DEG_results.P.Value':'pvalue',
            'DEG_results.adj.P.Val':'FDR',
            'AveExpr':'avg_expression',
            't':'tstat',
            'P.Value':'pvalue',
            'adj.P.Val':'FDR',
            'statistic':'tstat',
            'estimate':'log2FC',
            'pval':'pvalue',
            'std.error':'log2FC_se',
            'log2FoldChange':'log2FC',
            'lfcSE':'log2FC_se',
            'stat':'tstat',
            'padj':'FDR',
            'logFC':'log2FC'}

keep_cols = ['method', 'comparison', 'gene', 'log2FC', 'log2FC_se', 'pvalue', 'FDR', 'tstat', 'se', 'avg_expression', 'filename']

res = []
for m in methods.keys():
    for comp in files[m].keys():
        fn = opj(results_folder, methods[m], files[m][comp])
        tmp = pd.read_csv(fn)
        tmp = tmp.rename(col_map, axis=1)
        if not 'gene' in tmp.columns:
            tmp = tmp.rename({'Unnamed: 0':'gene'}, axis=1)
        tmp = tmp.assign(method=m, comparison=comp, filename=fn)
        res.append(tmp[[c for c in keep_cols if c in tmp.columns]])
res = pd.concat(res, axis=0)

res.to_csv(opj(results_folder, 'agregated_results_2023-MAR-15.csv'), index=False)

