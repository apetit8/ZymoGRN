#Infer networks/module from latidudinal data
#WGCNA setup####################################################################
library(WGCNA)
library(DESeq2)
library(tximport)
library(tximportData)
source("scripts/functions/RNAseq_analyses.R")
options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)
################################################################################
################################################################################
exp_counts <- read.csv("All_Data/genes_count_uncorrected.csv", sep=",")[,c(1:8,45:50)]
lat_counts2 <- read.csv("scripts/data/lat_gene_expr_raw_uncorrected.csv")
#
lat_counts2 <- lat_counts2[lat_counts2$X %in% exp_counts$X,]
exp_counts <- exp_counts[exp_counts$X %in% lat_counts2$X,]
#
data <- as.data.frame.matrix((cbind(lat_counts2[,2:31],round(2^(exp_counts[,2:ncol(exp_counts)])))))
row.names(data) <- exp_counts$X
################################################################################
#Import count data
dds <- DESeqDataSetFromMatrix(data, colData = as.data.frame(colnames(data)), design = ~ 0)

################################################################################
expression <-t(counts(estimateSizeFactors(dds), normalized=TRUE))
#As recommended, only keep genes that are strongly expressed :
#quant <- unlist(lapply(1:ncol(expression), FUN=function(i)( quantile(expression[,i], prob=c(.75), type=1) )))
#expression <- expression[,quant>15]
#expression <- log2(expression+1)
################################################################################
nGenes <- ncol(expression)
nSamples <- nrow(expression)
################################################################################
#2.a One-step automatic network construction and module detection
# Choose a set of soft-thresholding powers
powers <- c(c(1:15), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(expression, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
abline(h=0.9,col="red")

#Take appropriate power and apply it there 
#softPower <- 12
softPower <- 9
################################################################################
#Dealing with a big data set: Constructing a block-wise network and detecting modules
cor <- WGCNA::cor #Avoid a common bug
netwk <- blockwiseModules(expression, maxBlockSize = 5000, networkType = "unsigned",
                          power = softPower, TOMType = "signed", minModuleSize = 100,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE,
                          verbose = 3)
cor<-stats::cor
mergedColors <- labels2colors(netwk$colors)
moduleLabels <- netwk$colors
# Re-labeling blockwise modules
bwLabels <- matchLabels(mergedColors, moduleLabels)
#Converting the labels in colors to plot
bwmergedColors <- labels2colors(bwLabels)
probes <- colnames(expression)
moduleColors <- labels2colors(moduleLabels)
################################################################################
library(dplyr)
#Relating modules to characteristics and identifying important genes
module_df <- data.frame(gene_id = names(netwk$colors), colors = labels2colors(netwk$colors))
#Save in gene module the genes are
write.csv(module_df, "scripts/data/FULLDATA2_genes_x_module.csv")
#count_lat <- counts(estimateSizeFactors(dds), normalized=TRUE)
################################################################################