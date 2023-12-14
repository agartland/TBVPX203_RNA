
require(dplyr)
library(limma)
library(stringr)
library(magrittr)
library(edgeR)
library(limma)


data_filtering_EdgeR <- function(count_matrix)
{
  cpm <- cpm(count_matrix)
  lcpm <- cpm(count_matrix, log=TRUE)
  # fitering  genes
  keep <- rowSums(cpm(count_matrix)>1) >= 15
  d0 <- count_matrix[keep,]
  nsamples <- ncol(d0)
  
  
  return(d0)
}


data_normalization_edgeR <- function(filtered_count_matrix)
{
  # creating a DGE object
  normalized_matrix <- DGEList(filtered_count_matrix)
  normalized_matrix <- calcNormFactors(normalized_matrix,method = "TMM")
  return(normalized_matrix)
  
}



finding_DEGs <- function(meta_data,counts,time_list,treatment_group,model_parameter)
{
  sample_names = meta_data %>% 
    dplyr::filter(day %in% time_list) %>%
    dplyr::filter(Treatment_Group %in% treatment_group) %>%
    pull(samplename)
  
  
  mts_ = meta_data[sample_names,]
  cts_ = counts[,sample_names]
  
  
  # filtering and normalization of the dataset
  data1 <- data_filtering_EdgeR(cts_)
  d_ss <-  data_normalization_edgeR(data1)
  
  day <- mts_[["day"]]
  day <- as.character(day)
  model_matrix <-  model.matrix(~day, d_ss$samples)
  design <- model_matrix
  
  
  v <- voom(d_ss, design)
  
  
  
  v$targets = v$targets %>%
    mutate(day = as.factor(mts_[["day"]])) %>%
    mutate(sex = mts_[["sex"]]) %>%
    mutate(libID=colnames(cts_))%>%
    mutate(ptid=mts_[["ptid"]])
  
  
  klm     <- kmFit(dat = v,
                   model = model_parameter,
                   run.lme = TRUE,
                   use.weights = TRUE,
                   metrics = TRUE,
                   patientID = "ptid",
                   run.contrast = TRUE,
                   processors=3) 
  
  
  return(klm)
  
}




counts <- read.csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/salmon.merged.gene_counts_corrected.tsv",row.names = 1)
meta_data <- read.csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/salmon_metadata.csv", row.names = 1)
treatment_file <- read.csv("/fh/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/trt_pubid_2023-NOV-28.csv")

colnames(counts)<-gsub("X","P",colnames(counts))
meta_data$ptid <- str_replace(meta_data$ptid, "203", "")
meta_data$ptid <- str_replace(meta_data$ptid, "-", "")
meta_data$ptid <- str_replace(meta_data$ptid, "-", "_")


meta_data = merge(meta_data,treatment_file, by = "ptid")
meta_data$ptid = str_replace(meta_data$ptid, "_", "")
rownames(meta_data) <- meta_data$samplename



time_list <- list(c("0", "3"),c("0", "56"),c("56", "59"),c("56", "63"),c("56", "112"),c("56", "168"))
treatment_group <- c("2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)", "2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)")
model_parameter <- "~day+(1|ptid)+sex"




for (p in time_list) {
  nam <- paste("kim","_", paste0(p[1],"_",p[2]), sep = "")
  assign(nam, finding_DEGs(meta_data,counts,p,treatment_group,model_parameter))
}   




treatment_group <- c("2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)", 
                     "2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)",
                     "2 µg ID93 + 2 µg GLA-SE ", 
                     "2 µg ID93 + 2 µg GLA-SE")


kim_0_3_all4_group <- finding_DEGs(meta_data,counts,time_list = c("0", "3"),
                                   treatment_group,model_parameter)


sample_names = meta_data %>% 
  dplyr::filter(day %in% c("0","3")) %>%
  dplyr::filter(Treatment_Group %in% treatment_group) %>%
  pull(samplename)


mts_ = meta_data[sample_names,]
cts_ = counts[,sample_names]


# filtering and normalization of the dataset
data1 <- data_filtering_EdgeR(cts_)
d_ss <-  data_normalization_edgeR(data1)

day <- mts_[["day"]]
day <- as.character(day)
model_matrix <-  model.matrix(~day, d_ss$samples)
design <- model_matrix


v <- voom(d_ss, design)

model_parameter <- "~day+(1|ptid)+sex+adj"

v$targets = v$targets %>%
  mutate(day = as.factor(mts_[["day"]])) %>%
  mutate(sex = mts_[["sex"]]) %>%
  mutate(libID =colnames(cts_))%>%
  mutate(ptid = mts_[["ptid"]]) %>%
  mutate(adj = if_else(grepl("5", mts_[["Treatment_Group"]]),1,0))



kim_0_3_all4_group_adj <- kmFit(dat = v,
             model = model_parameter,
             run.lme = TRUE,
             use.weights = TRUE,
             metrics = TRUE,
             patientID = "ptid",
             run.contrast = TRUE,
             processors=3) 

kim_0_3_a <- kim_0_3$lme
kim_0_56_a  <- kim_0_56$lme
kim_56_59_a <- kim_56_59$lme
kim_56_63_a <- kim_56_63$lme
kim_56_112_a <- kim_56_112$lme
kim_56_168_a <-kim_56_168$lme
kim_0_3_all4_group_a <- kim_0_3_all4_group$lme
kim_0_3_all4_group_adj_a <- kim_0_3_all4_group_adj$lme



kim_0_3_b <- kim_0_3$lme.contrast
kim_0_56_b  <- kim_0_56$lme.contrast
kim_56_59_b <- kim_56_59$lme.contrast
kim_56_63_b <- kim_56_63$lme.contrast
kim_56_112_b <- kim_56_112$lme.contrast
kim_56_168_b <-kim_56_168$lme.contrast
kim_0_3_all4_group_b <- kim_0_3_all4_group$lme.contrast
kim_0_3_all4_group_adj_b <- kim_0_3_all4_group_adj$lme.contrast



kim_0_3_c <- kim_0_3$lme.fit
kim_0_56_c  <- kim_0_56$lme.fit
kim_56_59_c <- kim_56_59$lme.fit
kim_56_63_c <- kim_56_63$lme.fit
kim_56_112_c <- kim_56_112$lme.fit
kim_56_168_c <-kim_56_168$lme.fit
kim_0_3_all4_group_c <- kim_0_3_all4_group$lme.fit
kim_0_3_all4_group_adj_c <- kim_0_3_all4_group_adj$lme.fit

my_list <- list(kim_0_3_a,kim_0_56_a,kim_56_59_a,kim_56_63_a,kim_56_112_a,kim_56_168_a,
                kim_0_3_all4_group_a,kim_0_3_all4_group_adj_a)

# name the data frames

names(my_list) <- list("kim_0_3_lme","kim_0_56_lme","kim_56_59_lme","kim_56_63_lme",
                       "kim_56_112_lme","kim_56_168_lme","kim_0_3_all4_group_lme","kim_0_3_all4_group_adj_lme")

# save each new data frame as an individual .csv file based on its name

lapply(1:length(my_list), function(i) write.csv(my_list[[i]], 
                                                file = paste0(names(my_list[i]), ".csv"),
                                                row.names = FALSE))





my_list <- list(kim_0_3_b,kim_0_56_b,kim_56_59_b,kim_56_63_b,kim_56_112_b,kim_56_168_b,
                kim_0_3_all4_group_b,kim_0_3_all4_group_adj_b)

# name the data frames

names(my_list) <- list("kim_0_3_lme_contrast","kim_0_56_lme_contrast",
                       "kim_56_59_lme_contrast","kim_56_63_lme_contrast",
                       "kim_56_112_lme_contrast","kim_56_168_lme_contrast",
                       "kim_0_3_all4_group_lme_contrast","kim_0_3_all4_group_adj_lme_contrast")

# save each new data frame as an individual .csv file based on its name

lapply(1:length(my_list), function(i) write.csv(my_list[[i]], 
                                                file = paste0(names(my_list[i]), ".csv"),
                                                row.names = FALSE))




my_list <- list(kim_0_3_c,kim_0_56_c,kim_56_59_c,kim_56_63_c,kim_56_112_c,kim_56_168_c,
                kim_0_3_all4_group_c,kim_0_3_all4_group_adj_c)

# name the data frames

names(my_list) <- list("kim_0_3_lme_fit","kim_0_56_lme_fit",
                       "kim_56_59_lme_fit","kim_56_63_lme_fit",
                       "kim_56_112_lme_fit","kim_56_168_lme_fit",
                       "kim_0_3_all4_group_lme_fit","kim_0_3_all4_group_adj_lme_fit")

# save each new data frame as an individual .csv file based on its name

lapply(1:length(my_list), function(i) write.csv(my_list[[i]], 
                                                file = paste0(names(my_list[i]), ".csv"),
                                                row.names = FALSE))
