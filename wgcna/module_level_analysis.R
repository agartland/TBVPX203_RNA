

library(kimma)
require(dplyr)
library(limma)
library(stringr)
library(lme4)
library(tidyverse)
library(lmerTest)
library(broom.mixed)





meta_data <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_metadata.csv", row.names = 1)
treatment_file <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/trt_pubid_2023-NOV-28.csv")


meta_data$ptid <- str_replace(meta_data$ptid, "203", "")
meta_data$ptid <- str_replace(meta_data$ptid, "-", "")
meta_data$ptid <- str_replace(meta_data$ptid, "-", "_")


data_modeling <- function(f,master_model_matrix_, formula,time,model_name)
{
  test_data <- master_model_matrix_  %>%
    dplyr::filter(day %in%  time) %>%
    dplyr::filter(module == f)
  
  mod4 <-lmer(form4, data = test_data)
  mod4_r <- broom.mixed::tidy(mod4,conf.int = TRUE)  %>%
    dplyr::filter(effect == "fixed") %>%
    dplyr::select(-c("effect","group"))  %>% 
    dplyr::filter(term != "(Intercept)")
  mod4_AIC <- lmer(form4, data = test_data,  REML = TRUE)
  mod4_r_AIC <- glance(mod4_AIC)%>%
    dplyr::select(c("AIC","BIC"))
  
  results <- cbind(mod4_r,mod4_r_AIC)
  results["module"]<-f
  results["model"]<-model_name
  
  return(results)
}
meta_data = merge(meta_data,treatment_file, by = "ptid")
meta_data$ptid = str_replace(meta_data$ptid, "_", "")
rownames(meta_data) <- meta_data$samplename
treatment_group <- c("2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)", 
                     "2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)")


sample_A = meta_data %>% 
  dplyr::filter(Treatment_Group %in% treatment_group) %>%
  pull(samplename)


wgcna_file <- read.csv("X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/21Nov2023/module_results/module_norm_cts.csv",row.names = 1)
wgcna_file = wgcna_file %>% 
  dplyr::filter(sampleid %in% sample_A)




master_model_matrix_ = wgcna_file %>%
  mutate(day_59 = if_else(str_sub(sampleid, - 2, - 1) == "59", 1, 0)) %>%
  mutate(day_63 = if_else(str_sub(sampleid, - 2, - 1) == "63", 1, 0)) %>%
  mutate(day_56 = if_else(str_sub(sampleid, - 2, - 1) == "56", 1, 0))%>%
  mutate(dose_2 = if_else(str_sub(Treatment_Group,-21,-21) == "2", 1, 0)) %>%
  mutate(dose_3 = if_else(str_sub(Treatment_Group,-21,-21) == "3", 1, 0)) %>%
  mutate(sex__ = if_else(sex_ == "female", 1, 0))


module = list("green","blue","turquoise","red","black","yellow","brown","pink")
model_name <- "no_lmer_interactions"
form4 = as.formula(paste("eigengene", "~", "day_59","+","dose_2","+","sex_","+", "(1|ptid)"))
time <- c("59","56")



store = list()
store_sum = list()
for (f in module)
{
  pval4 = list()
  x_1 <- data_modeling(f,master_model_matrix_, form4,time,model_name)
  tag1 = paste0(f,"no_lmer_interactions")
  store[[tag1]]<-x_1

}

final_dta_4<-bind_rows(store)
day_59 <-rbind(final_dta_1,final_dta_2,final_dta_3,final_dta_4)


write.csv(day_59,"day_59_module_level_analysis.csv")








