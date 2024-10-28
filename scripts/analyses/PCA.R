source("scripts/functions/RNAseq_analyses.R")
################################################################################
library(ggplot2)
library(ggfortify)
################################################################################
uncorrected_counts <- read.csv("All_Data/genes_count_uncorrected.csv", sep=",")[,2:80]
corrected_counts <- read.csv("All_Data/genes_count_corrected.csv", sep=",")[,2:80]

df.PCA <- as.data.frame.matrix(t(uncorrected_counts))
df.PCA$Pop <- rownames(df.PCA)
df.PCA$Temp <- stringr::str_split(rownames(df.PCA), pattern="\\.",n=5 , simplify = TRUE)[,3]
df.PCA$Strain <- stringr::str_split(rownames(df.PCA), pattern="\\.",n=5 , simplify = TRUE)[,1]
df.PCA$Strain <- stringr::str_replace_all(df.PCA$Strain, "S01", "MGGP01")
df.PCA$Strain <- stringr::str_replace_all(df.PCA$Strain, "S44", "MGGP44")
df.PCA$Year <- stringr::str_split(rownames(df.PCA), pattern="\\.",n=5 , simplify = TRUE)[,5]

pca_rest1 <- prcomp(df.PCA[,1:(ncol(df.PCA)-4)], scale. = FALSE)

gg1 <- autoplot(pca_rest1, data = df.PCA, colour = 'Year', shape='Strain', size = 3, label = FALSE, x = 1, y = 2,
                loadings = FALSE, cex=2, loadings.label = FALSE, loadings.label.size = 4)+  
  ggtitle("PCA on gene expression (batch uncorrected)")+theme_bw() +
  scale_color_manual(values =  c("khaki","orange"))

gg2 <- autoplot(pca_rest1, data = df.PCA, colour = 'Temp', shape='Strain', size = 3, label = FALSE, x = 3, y = 4,
                loadings = FALSE, cex=.8, loadings.label = FALSE, loadings.label.size = 4, loadings.label.colour = "darkgrey", loadings.colour = "lightgrey")+  
  ggtitle("PCA on gene expression (batch uncorrected)")+theme_bw()  +
  scale_color_manual(values =  c("lightskyblue","lightcoral"))

###
corrected_counts <- read.csv("All_Data/genes_count_corrected.csv", sep=",")[,2:80]
rownames(corrected_counts) <- read.csv("All_Data/genes_count_corrected.csv", sep=",")[,1]

df.PCA2 <- as.data.frame.matrix(t(corrected_counts))
df.PCA2$Pop <- rownames(df.PCA2)
#df.PCA2$Batch <- df.PCA$Batch
df.PCA2$Temp <- stringr::str_split(rownames(df.PCA2), pattern="\\.",n=5 , simplify = TRUE)[,3]
df.PCA2$Strain <- stringr::str_split(rownames(df.PCA2), pattern="\\.",n=5 , simplify = TRUE)[,1]
df.PCA2$Strain <- stringr::str_replace_all(df.PCA2$Strain, "S01", "MGGP01")
df.PCA2$Strain <- stringr::str_replace_all(df.PCA2$Strain, "S44", "MGGP44")
df.PCA2$Year <- stringr::str_split(rownames(df.PCA), pattern="\\.",n=5 , simplify = TRUE)[,5]
df.PCA2$Lines <- stringr::str_split(rownames(df.PCA2), pattern="\\.", simplify = TRUE)[,2] #paste0(stringr::str_split(rownames(df.PCA2), pattern="\\.", simplify = TRUE)[,1],stringr::str_split(rownames(df.PCA2), pattern="\\.", simplify = TRUE)[,2])
df.PCA2$Lines <- stringr::str_replace_all(df.PCA2$Lines, "AAA", "Ancestor")

pca_rest2 <- prcomp(df.PCA2[,1:(ncol(df.PCA2)-5)], scale. = FALSE)

gg3 <- autoplot(pca_rest2, data = df.PCA2, colour = 'Temp', shape='Strain', size = 3, label = FALSE, x = 1, y = 2, 
                loadings = FALSE, alpha=1, cex=.8, loadings.label = FALSE, loadings.label.size = 4, loadings.label.colour = "darkgrey", loadings.colour = "lightgrey")+  
  ggtitle("PCA on gene expression (batch corrected)")+theme_bw() +
  scale_color_manual(values =  c("lightskyblue","lightcoral"))

gg4 <- autoplot(pca_rest2, data = df.PCA2, colour = 'Temp', shape='Strain', size = 3, label = FALSE, x = 3, y = 4, 
                loadings = FALSE, cex=.8, loadings.label = FALSE, loadings.label.size = 4, loadings.label.colour = "darkgrey", loadings.colour = "lightgrey")+  
  ggtitle("PCA on gene expression (batch corrected)")+theme_bw()  +
  scale_color_manual(values =  c("lightskyblue","lightcoral"))

gg5 <- autoplot(pca_rest2, data = df.PCA2, colour = 'Lines', shape='Strain', size = 3, label = FALSE, x = 3, y = 4,
                loadings = FALSE, cex=0.8, loadings.label = FALSE, loadings.label.size = 4, loadings.label.colour = "darkgrey", loadings.colour = "lightgrey")+  
  ggtitle("PCA on gene expression (batch corrected)")+theme_bw() +
  scale_color_manual(values =  c("royalblue1","royalblue","blue3","orangered","red","red3","grey20","pink","pink2","pink3",
                                 "royalblue1","royalblue","orangered","red","red3","grey20","pink","pink2","pink3")) #c("royalblue1","royalblue","blue3","orangered","red","red3","grey50","mediumpurple","purple","darkorchid","skyblue1","skyblue2","orange","sandybrown","orange3","grey","pink","pink2","pink3")


pdf(paste0("figures/PCA_batch_correction.pdf"), width = 5.5, height = 4)
gg1
gg2
gg3
gg4
gg5
dev.off()

