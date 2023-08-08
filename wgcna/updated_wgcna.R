library(WGCNA)
library(dplyr)
library(DEG)
library(stringr)
library(limma)
library(edgeR)
library(ggplot2)
library(MASS)
library(Hmisc)



# normalized counts for kimma DEG output with the associated samples

dat <- readRDS("WGCNA_stricter_DEG_expression_dat.rds")
datExpr<-as.data.frame(dat)
datExpr0 = as.data.frame(t(datExpr))


####################################### WGCNA ##############################################

sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2)



# Plot a line to show the cut
abline(h = 150, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)




# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.850,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# signed correlation for expression data #
softPower = 1
corr <- rcorr(as.matrix(datExpr))$r
adj_2<- ((corr+1)/2)^5
adjacency = adjacency(datExpr, power = 1)


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adj_2);
dissTOM = 1-TOM



# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);



# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 10;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")




# Rename to moduleColors
moduleColors = dynamicColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;


load("WGCNA_modules_DEG_expression_dat.RData")


a=list()

color_s <- unique( moduleColors)
for (i in color_s){
  a[[i]]<-names(datExpr)[ moduleColors == i]
}
b<-c()
for( i in 1:length(a))
{
  b[i]=max(length(a[[i]]))
}
b=max(b)


for( i in 1:length(a))
{
  a[[i]]<-c(a[[i]],rep(NA, b - length(a[[i]])))
}



# Append columns within for loop
data <- data.frame(a)
names(data) <- color_s











####################################### PRINCIPLE COMPONENT ANALYSIS ##############################################


# load the normalized dataset across all tyhe sample 
counts_file = 'R_cts.csv'
meta_file = 'R_meta.csv'
cts = read.csv(counts_file, row.names = 1)
mts = read.csv(meta_file)
rownames(mts) = mts[['pubvid']]
# select appropriate samples
y <- DGEList(cts)
#/ calculate TMM normalization factors:
y <- calcNormFactors(y)
#/ get the normalized counts:
cpms <- cpm(y, log=FALSE)


#find the standard deviation across all the samples using normalized data across all the samples 
data1<-cpms
data_stan<- t(data1)
data_sd2<-apply(data_stan,2,function(x) (x - mean(x)) / sd(x))
data_sd2<-data_sd2[ , colSums(is.na(data_sd2))==0]





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



modules <- module_data
module<-modules$turquoise
second_module <- na.omit(modules$brown)
third_module <- na.omit(modules$blue)


turquoise <- data_pca(cpms_bcg,data_sd2,module)
brown <- data_pca(cpms_bcg,data_sd2,second_module)
blue <- data_pca(cpms_bcg,data_sd2, third_module)


eigen_genes <- cbind(turquoise,brown,blue)
write.csv(eigen_genes,"/fh/fast/gilbert_p/hsrivast/external_test/TB602_results/wgcna_eigen_gene_results/eigen_genes_bcg.csv")

