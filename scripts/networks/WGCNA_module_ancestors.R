#Module inference with WGCNA
#WGCNA setup####################################################################
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)
################################################################################
################################################################################
corrected_counts <- read.csv("All_Data/genes_count_corrected.csv", sep=",")
#
expression <- t(corrected_counts[,c(2:8, 45:50)])
colnames(expression) <- corrected_counts[,1]
#
quant <- unlist(lapply(1:ncol(expression), FUN=function(i)( quantile(expression[,i], prob=c(.75), type=1) )))
#
nGenes <- ncol(expression)
nSamples <- nrow(expression)
################################################################################
#2.a One-step automatic network construction and module detection
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
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
softPower <- 11
################################################################################
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
write.csv(module_df, "scripts/data/Ancestors_genes_x_module_corrected_exp_data.csv")