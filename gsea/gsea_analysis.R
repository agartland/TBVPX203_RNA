install.packages("remotes")
remotes::install_github("saezlab/COSMOS")
install.packages("msigdbr")
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(tidyr)
library(COSMOS)

kim_0_3 <- read.csv("kim_0_3.csv")
kim_0_56 <- read.csv("kim_0_56.csv")
kim_56_59 <- read.csv("kim_56_59.csv")
kim_56_63 <- read.csv("kim_56_63.csv")
kim_56_112 <- read.csv("kim_56_112.csv")
kim_56_168 <- read.csv("kim_56_168.csv")
kim_0_3_all4_group <- read.csv("kim_0_3_all4_group.csv")
kim_0_3_all4_group_adj <- read.csv("kim_0_3_all4_group_adj.csv")



my_list_genes <- unique(c(kim_0_3$gene,kim_0_56$gene,kim_56_59$gene,kim_56_63$gene,kim_56_112$gene,                          
                          kim_56_168$gene,kim_0_3_all4_group$gene,kim_0_3_all4_group_adj$gene))
all_genes <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/salmon.merged.gene_counts_corrected.tsv")
genes <- all_genes$gene_id

wgcna <- read.csv("new_wgcna_module.csv")

deg_gene_list_1 <- na.omit(wgcna$turquoise)
deg_gene_list_2 <- na.omit(wgcna$brown)
deg_gene_list_3 <- na.omit(wgcna$yellow)
deg_gene_list_4 <- na.omit(wgcna$green)
deg_gene_list_5 <- na.omit(wgcna$blue)
deg_gene_list_6 <- na.omit(wgcna$red)
deg_gene_list_7 <- na.omit(wgcna$black)
deg_gene_list_8 <- na.omit(wgcna$pink)


m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, human_gene_symbol)
head(m_t2g)

gmt.file <- file.path("X:/fast/gilbert_p/fg_data/BTM/BTM_for_GSEA.gmt")
pathways <- gmtPathways(gmt.file)


datalist = list()
pathwayData <- data.frame(pathway= character(0), gene= character(0))
for(i in 1:length(pathways))
{
  df <- data.frame(matrix(unlist(pathways[[i]]), nrow=length(pathways[[i]]), byrow=TRUE))
  colnames(df) <- "human_gene_symbol"
  df['gs_names'] <- names(pathways[i])
  datalist[[i]] <- df 
  
}
pathwayData <- do.call(rbind, datalist)
m2g_gmt <- pathwayData["gs_names"]
m2g_gmt["human_gene_symbol"] <- pathwayData["human_gene_symbol"]





em_C7_deg_vs_all <- enricher(my_list_genes, universe=genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_turquoise_vs_deg <- enricher(deg_gene_list_1, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_brown_vs_deg <- enricher(deg_gene_list_2, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_yellow_vs_deg <- enricher(deg_gene_list_3, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_green_vs_deg <- enricher(deg_gene_list_4, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_blue_vs_deg <- enricher(deg_gene_list_5, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_red_vs_deg <- enricher(deg_gene_list_6, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_black_vs_deg <- enricher(deg_gene_list_7, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)
em_C7_pink_vs_deg <- enricher(deg_gene_list_8, universe=my_list_genes, minGSSize=5, TERM2GENE=m_t2g)






em_BTM_deg_vs_all <- enricher(my_list_genes, universe=genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_turquoise_vs_deg <- enricher(deg_gene_list_1, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_brown_vs_deg <- enricher(deg_gene_list_2, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_yellow_vs_deg <- enricher(deg_gene_list_3, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_green_vs_deg <- enricher(deg_gene_list_4, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_blue_vs_deg <- enricher(deg_gene_list_5, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_red_vs_deg <- enricher(deg_gene_list_6, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_black_vs_deg <- enricher(deg_gene_list_7, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)
em_BTM_pink_vs_deg <- enricher(deg_gene_list_8, universe=my_list_genes, minGSSize=5, TERM2GENE=m2g_gmt)


my_list <- list(em_C7_deg_vs_all,em_C7_turquoise_vs_deg,em_C7_brown_vs_deg,
                em_C7_yellow_vs_deg,em_C7_green_vs_deg,
                em_C7_blue_vs_deg,em_C7_red_vs_deg,em_C7_black_vs_deg,em_C7_pink_vs_deg,em_BTM_deg_vs_all,
                em_BTM_turquoise_vs_deg,em_BTM_brown_vs_deg,em_BTM_yellow_vs_deg,em_BTM_green_vs_deg,
                em_BTM_blue_vs_deg,em_BTM_red_vs_deg,em_BTM_black_vs_deg,em_BTM_pink_vs_deg)

# name the data frames

names(my_list) <- list("em_C7_deg_vs_all","em_C7_turquoise_vs_deg","em_C7_brown_vs_deg","em_C7_yellow_vs_deg","em_C7_green_vs_deg",
                       "em_C7_blue_vs_deg","em_C7_red_vs_deg","em_C7_black_vs_deg","em_C7_pink_vs_deg","em_BTM_deg_vs_all",
                       "em_BTM_turquoise_vs_deg","em_BTM_brown_vs_deg","em_BTM_yellow_vs_deg","em_BTM_green_vs_deg",
                       "em_BTM_blue_vs_deg","em_BTM_red_vs_deg","em_BTM_black_vs_deg","em_BTM_pink_vs_deg")

# save each new data frame as an individual .csv file based on its name

lapply(1:length(my_list), function(i) write.csv(my_list[[i]], 
                                                file = paste0(names(my_list[i]), ".csv"),
                                                row.names = FALSE))



