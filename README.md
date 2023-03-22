# TBVPX-203 transcriptomic analysis

## Code folders

 - `deg`: scripts for differential gene expression (will depend on `himangi2003/friendly_waddle`)
 - `gsea`: scripts for gene set enrichment analysis
 - `wgcna`: scripts for weighted gene correlation network analysis
 - `ircorr`: scripts to evaluate correlation of immune responses with gene modules


## Data folders

Folder containing raw immune response and clinical data:
`T:\vaccine\IDRI\TBVPX203\qdata`
`qdata_folder = opj(_trials, 'vaccine', 'IDRI', 'TBVPX203', 'qdata')`

Folder containing DEG results from all methods:
`X:\fast\gilbert_p\fg_data\SCRI\TBVPX-203\RNA\2019Dec\Results`

Salmon counts file:
`X:\fast\gilbert_p\fg_data\SCRI\TBVPX-203\RNA\salmon_merged_gene_counts.csv`

Treatment assignments file:
`X:\fast\gilbert_p\fg_data\SCRI\TBVPX-203\RNA\trt_pubid_2022-DEC-19.csv`

BTM file for GSEA:
`X:\fast\gilbert_p\fg_data\BTM\BTM_for_GSEA_20131008.gmt`