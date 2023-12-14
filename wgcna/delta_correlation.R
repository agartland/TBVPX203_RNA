
library(kimma)
require(dplyr)
library(limma)
library(stringr)
library(lme4)
library(tidyverse)
library(lmerTest)

#data loading and preprocessing
counts <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/salmon.merged.gene_counts_corrected.tsv",row.names = 1)
meta_data <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_metadata.csv", row.names = 1)
treatment_file <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/trt_pubid_2023-NOV-28.csv")
norm_counts <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/log_normalized_counts.csv",row.names = 1)
colnames(norm_counts )<-gsub("X","P",colnames(norm_counts ))


colnames(counts)<-gsub("X","P",colnames(counts))
meta_data$ptid <- str_replace(meta_data$ptid, "203", "")
meta_data$ptid <- str_replace(meta_data$ptid, "-", "")
meta_data$ptid <- str_replace(meta_data$ptid, "-", "_")


meta_data = merge(meta_data,treatment_file, by = "ptid")
meta_data$ptid = str_replace(meta_data$ptid, "_", "")
rownames(meta_data) <- meta_data$samplename
treatment_group <- c("2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)", 
                     "2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)")


kim_0_3 <- read.csv("degs/kim_0_3.csv")
kim_0_56 <- read.csv("degs/kim_0_56.csv")
kim_56_59 <- read.csv("degs/kim_56_59.csv")
kim_56_63 <- read.csv("degs/kim_56_63.csv")
kim_56_112 <- read.csv("degs/kim_56_112.csv")
kim_56_168 <- read.csv("degs/kim_56_168.csv")
kim_0_3_all4_group <- read.csv("degs/kim_0_3_all4_group.csv")
kim_0_3_all4_group_adj <- read.csv("degs/kim_0_3_all4_group_adj.csv")



my_list_genes <- unique(c(kim_0_3$gene,kim_0_56$gene,kim_56_59$gene,kim_56_63$gene,kim_56_112$gene,                          
                          kim_56_168$gene,kim_0_3_all4_group$gene,kim_0_3_all4_group_adj$gene))





norm_counts_deg <-norm_counts[my_list_genes,] 

treatment_group <- c("2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)", 
                     "2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)")


wgcna <- read.csv("new_wgcna_module.csv")

deg_gene_list_1 <- na.omit(wgcna$turquoise)
deg_gene_list_2 <- na.omit(wgcna$brown)
deg_gene_list_3 <- na.omit(wgcna$yellow)
deg_gene_list_4 <- na.omit(wgcna$green)
deg_gene_list_5 <- na.omit(wgcna$blue)
deg_gene_list_6 <- na.omit(wgcna$red)
deg_gene_list_7 <- na.omit(wgcna$black)
deg_gene_list_8 <- na.omit(wgcna$pink)


################################################## 0 3 difference ##########################################################################################
time_list <- "0"
ptid_0 = meta_data %>%
  dplyr::filter(day %in% time_list) %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


time_list <- "3"
ptid_3 = meta_data %>%
  dplyr::filter(day %in% time_list) %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_0_3 <- as.character(ptid_0[(ptid_0 %in% ptid_3)])
ptid_0_3_ni <- as.character(ptid_0[!(ptid_0 %in% ptid_3)])


dt_n <-data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_0_3_ni)))
rownames(dt_n) <- my_list_genes
colnames(dt_n) <- ptid_0_3_ni
dt <- data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_0_3)))
rownames(dt) <- my_list_genes
colnames(dt) <- ptid_0_3
for ( i in ptid_0_3)
{
  sample_name_0 = meta_data %>%
    dplyr::filter(ptid == i) %>%
    dplyr::filter(day %in% c("0")) %>%
    pull(samplename)
  
  sample_name_3 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("3")) %>%
    pull(samplename)
  
  time_3_sample <- norm_counts_deg[,sample_name_3]
  time_0_sample <- norm_counts_deg[,sample_name_0]
  dt[i]<-c(time_3_sample-time_0_sample)

}

dt_0_3<-cbind(dt,dt_n)

################################################## 56 59 difference ##########################################################################################

ptid_56 = meta_data %>%
  dplyr::filter(day %in% "56") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_59 = meta_data %>%
  dplyr::filter(day %in% "59") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_56_59 <- ptid_59[(ptid_59 %in% ptid_56)]
ptid_56_59_ni <- ptid_59[!(ptid_59 %in% ptid_56)]


dt_n <-data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_59_ni)))
rownames(dt_n) <- my_list_genes
colnames(dt_n) <- ptid_56_59_ni
dt <- data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_59)))
rownames(dt) <- my_list_genes
colnames(dt) <- ptid_56_59
for ( i in ptid_56_59)
{
  sample_name_56 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("56")) %>%
    pull(samplename)
  
  sample_name_59 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("59")) %>%
    pull(samplename)
  
  time_59_sample <- norm_counts_deg[,sample_name_59]
  time_56_sample <- norm_counts_deg[,sample_name_56]
  dt[i]<-c(time_59_sample-time_56_sample)
  
}

dt_56_59<-cbind(dt,dt_n)

################################################## 56 59 difference ##########################################################################################

ptid_56 = meta_data %>%
  dplyr::filter(day %in% "56") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_63 = meta_data %>%
  dplyr::filter(day %in% "63") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_56_63 <- ptid_63[(ptid_63 %in% ptid_56)]
ptid_56_63_ni <- ptid_63[!(ptid_63 %in% ptid_56)]


dt_n <-data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_63_ni)))
rownames(dt_n) <- my_list_genes
colnames(dt_n) <- ptid_56_63_ni
dt <- data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_63)))
rownames(dt) <- my_list_genes
colnames(dt) <- ptid_56_63
for ( i in ptid_56_63)
{
  sample_name_56 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("56")) %>%
    pull(samplename)
  
  sample_name_63 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("63")) %>%
    pull(samplename)
  
  time_63_sample <- norm_counts_deg[,sample_name_63]
  time_56_sample <- norm_counts_deg[,sample_name_56]
  dt[i]<-c(time_63_sample-time_56_sample)
  
}

dt_56_63<-cbind(dt,dt_n)

###############################################################################################################################################################################################################


module <- c("turquoise","brown","yellow","green","blue","red","black","pink")
deg_modules <- list(deg_gene_list_1,deg_gene_list_2,deg_gene_list_3,
                    deg_gene_list_4,deg_gene_list_5,deg_gene_list_6, 
                    deg_gene_list_7,deg_gene_list_8)

data_n <-list()
for(i in 1:length(deg_modules))
  
{ 
  gene = na.omit(na.omit(deg_modules[[i]]))
  dt_0_3_mod <- dt_0_3[gene,]
  mod_0_3<-colMeans(dt_0_3_mod)
  data_n[[i]] <-c(mod_0_3) 
  
}

delta_A_0_3 <- as.data.frame(do.call(rbind, data_n))
rownames(delta_A_0_3)<-module


data_n <-list()
for(i in 1:length(deg_modules))
  
{ 
  gene = na.omit(na.omit(deg_modules[[i]]))
  dt_56_59_mod <- dt_0_3[gene,]
  mod_56_59<-colMeans(dt_56_59_mod)
  data_n[[i]] <-c(mod_56_59) 
  
}

delta_A_56_59 <- as.data.frame(do.call(rbind, data_n))
rownames(delta_A_56_59)<-module



data_n <-list()
for(i in 1:length(deg_modules))
  
{ 
  gene = na.omit(na.omit(deg_modules[[i]]))
  dt_56_63_mod <- dt_0_3[gene,]
  mod_56_63<-colMeans(dt_56_63_mod)
  data_n[[i]] <-c(mod_56_63) 
  
}

delta_A_56_63 <- as.data.frame(do.call(rbind, data_n))
rownames(delta_A_56_63)<-module

############################################################################################################################################################################################

################################################## 0 3 difference ##########################################################################################
time_list <- "0"
ptid_0 = meta_data %>%
  dplyr::filter(day %in% time_list) %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


time_list <- "3"
ptid_3 = meta_data %>%
  dplyr::filter(day %in% time_list) %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_0_3 <- as.character(ptid_0[(ptid_0 %in% ptid_3)])
ptid_0_3_ni <- as.character(ptid_0[!(ptid_0 %in% ptid_3)])


dt_n <-data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_0_3_ni)))
rownames(dt_n) <- my_list_genes
colnames(dt_n) <- ptid_0_3_ni
dt <- data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_0_3)))
rownames(dt) <- my_list_genes
colnames(dt) <- ptid_0_3
for ( i in ptid_0_3)
{
  sample_name_0 = meta_data %>%
    dplyr::filter(ptid == i) %>%
    dplyr::filter(day %in% c("0")) %>%
    pull(samplename)
  
  sample_name_3 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("3")) %>%
    pull(samplename)
  
  time_3_sample <- norm_counts_deg[,sample_name_3]
  time_0_sample <- norm_counts_deg[,sample_name_0]
  dt[i]<-c(time_3_sample-time_0_sample)

}

dt_0_3<-cbind(dt,dt_n)

################################################## 56 59 difference ##########################################################################################

ptid_56 = meta_data %>%
  dplyr::filter(day %in% "56") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_59 = meta_data %>%
  dplyr::filter(day %in% "59") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_56_59 <- ptid_59[(ptid_59 %in% ptid_56)]
ptid_56_59_ni <- ptid_59[!(ptid_59 %in% ptid_56)]


dt_n <-data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_59_ni)))
rownames(dt_n) <- my_list_genes
colnames(dt_n) <- ptid_56_59_ni
dt <- data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_59)))
rownames(dt) <- my_list_genes
colnames(dt) <- ptid_56_59
for ( i in ptid_56_59)
{
  sample_name_56 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("56")) %>%
    pull(samplename)
  
  sample_name_59 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("59")) %>%
    pull(samplename)
  
  time_59_sample <- norm_counts_deg[,sample_name_59]
  time_56_sample <- norm_counts_deg[,sample_name_56]
  dt[i]<-c(time_59_sample-time_56_sample)
  
}

dt_56_59<-cbind(dt,dt_n)

################################################## 56 59 difference ##########################################################################################

ptid_56 = meta_data %>%
  dplyr::filter(day %in% "56") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_63 = meta_data %>%
  dplyr::filter(day %in% "63") %>%
  dplyr::filter(Treatment_Group %in% treatment_group ) %>%
  pull(ptid)


ptid_56_63 <- ptid_63[(ptid_63 %in% ptid_56)]
ptid_56_63_ni <- ptid_63[!(ptid_63 %in% ptid_56)]


dt_n <-data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_63_ni)))
rownames(dt_n) <- my_list_genes
colnames(dt_n) <- ptid_56_63_ni
dt <- data.frame(matrix(nrow = length(my_list_genes), ncol = length(ptid_56_63)))
rownames(dt) <- my_list_genes
colnames(dt) <- ptid_56_63
for ( i in ptid_56_63)
{
  sample_name_56 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("56")) %>%
    pull(samplename)
  
  sample_name_63 = meta_data %>%
    dplyr::filter(ptid ==  i) %>%
    dplyr::filter(day %in% c("63")) %>%
    pull(samplename)
  
  time_63_sample <- norm_counts_deg[,sample_name_63]
  time_56_sample <- norm_counts_deg[,sample_name_56]
  dt[i]<-c(time_63_sample-time_56_sample)
  
}

dt_56_63<-cbind(dt,dt_n)

###############################################################################################################################################################################################################


module <- c("turquoise","brown","yellow","green","blue","red","black","pink")
deg_modules <- list(deg_gene_list_1,deg_gene_list_2,deg_gene_list_3,
                    deg_gene_list_4,deg_gene_list_5,deg_gene_list_6, 
                    deg_gene_list_7,deg_gene_list_8)

data_n <-list()
for(i in 1:length(deg_modules))
  
{ 
  gene = na.omit(na.omit(deg_modules[[i]]))
  dt_0_3_mod <- dt_0_3[gene,]
  mod_0_3<-colMeans(dt_0_3_mod)
  data_n[[i]] <-c(mod_0_3) 
  
}

delta_A_0_3 <- as.data.frame(do.call(rbind, data_n))
rownames(delta_A_0_3)<-module


data_n <-list()
for(i in 1:length(deg_modules))
  
{ 
  gene = na.omit(na.omit(deg_modules[[i]]))
  dt_56_59_mod <- dt_0_3[gene,]
  mod_56_59<-colMeans(dt_56_59_mod)
  data_n[[i]] <-c(mod_56_59) 
  
}

delta_A_56_59 <- as.data.frame(do.call(rbind, data_n))
rownames(delta_A_56_59)<-module



data_n <-list()
for(i in 1:length(deg_modules))
  
{ 
  gene = na.omit(na.omit(deg_modules[[i]]))
  dt_56_63_mod <- dt_0_3[gene,]
  mod_56_63<-colMeans(dt_56_63_mod)
  data_n[[i]] <-c(mod_56_63) 
  
}

delta_A_56_63 <- as.data.frame(do.call(rbind, data_n))
rownames(delta_A_56_63)<-module