####################################### PRINCIPLE COMPONENT ANALYSIS ##############################################

library(dplyr)

# load the normalized data counts for all the genes and all the samples
data1<-read.csv("normalized_data_all.csv",row.names = 1)
data_stan<- t(data1)
data_sd2<-apply(data_stan,2,function(x) (x - mean(x)) / sd(x))




#' function to find the eigen genes for a particular module 
#'
#' @param dat : normalized counts across the specific samples (subset of the 
#' samples used to find DEG GENES) and DEG genes
#' @param data_sd data obtained after calculating standard deviation 
#' across all the samples using normalized data across all the samples
#' @param first_module module for which the eigengenes need to be calculated
#'
#' @return eigengenes across all the samples for the module
data_pca <- function(dat,data_sd,module)
{
  data2_all_samples_first_module <- dplyr::select(as.data.frame(data_sd),
                                                  all_of(module))
  v1<-colnames(dat)
  data2 <- data2_all_samples_first_module[rownames
                                          (data2_all_samples_first_module)
                                          %in% v1, ]
  
  pca = prcomp(data2, scale = F, center = F)
  eg1 <- as.matrix(data2_all_samples_first_module) %*% pca$rotation
  eg1 <- eg1[,1]
  
  return(eg1)
}

# load the modules data
modules <-read.csv("WGCNA_modules_stricter_DEG_expression.csv")
first_module<-modules$turquoise
second_module <- na.omit(modules$brown)
third_module <- na.omit(modules$green)
fourth_module <- na.omit(modules$yellow)
fifth_module <- na.omit(modules$blue)
sixth_module <- na.omit(modules$red)
seventh_module <- na.omit(modules$black)

#calculate the first principole component for each module
turquoise <- data_pca(dat,data_sd2,first_module)
brown <- data_pca(dat,data_sd2,second_module)
green <- data_pca(dat,data_sd2,third_module)
yellow <- data_pca(dat,data_sd2, fourth_module)
blue <- data_pca(dat,data_sd2, fifth_module)
red <- data_pca(dat,data_sd2,sixth_module)
black <- data_pca(dat,data_sd2,seventh_module)



eigen_genes <- cbind(turquoise,brown,green,yellow,blue,red,black)
write.csv(eigen_genes,"eigen_genes.csv")

require(ggplot2)
ggplot(as.data.frame(pca$x), aes(x= PC1 , y= PC2 ) )+
  geom_point()

write.csv(network,"network_dta_stricter_DEG.csv")
plot(pca)
