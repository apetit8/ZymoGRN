#Infer networks/module from latidudinal data
#WGCNA setup####################################################################
library(WGCNA)
library(DESeq2)
library(tximport)
library(tximportData)
library(stringr)
source("scripts/functions/RNAseq_analyses.R")
options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)
################################################################################
################################################################################
#List of ribosomes from annotation file
ribosomes <- str_split(subset.pattern(read.csv(file="All_Data/Annotation_r.csv", sep='\t'), "ribosom", 4)[,3], as.character("\\."), n=2, simplify = TRUE)[,1]
min <-  10 #min = reads mean min for genes to be considered as outliers
#"Group" or coloumn of ref for ComBat_seq
group <- 8
################################################################################
#ONLY 01
################################################################################
#Data
quantData <- list.files(path = "All_Data/Latitudinal/", full.names=TRUE, recursive = TRUE, pattern = "t_data.ctab")
quantData <- quantData[c(1:12,19:36)]
quantData <- quantData[order(as.character(quantData), method = c("radix"))]
names(quantData) <- str_split(quantData, "/", simplify = TRUE)[,4] #select column corresponding to sample name
#Design info
samp.info <- read.csv("All_Data/Latitudinal/samples_design.csv", sep = ",")
rownames(samp.info) <- samp.info$X
samp.info <- samp.info[c(1:12,19:36),2:4]
samp.info[,1:ncol(samp.info)] <- lapply(samp.info[,1:ncol(samp.info)] , factor)
################################################################################
#Import count data
tmp <- read_tsv(quantData[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi1 <- tximport(files=quantData, type = "stringtie", tx2gene = tx2gene, readLength=75) 
dds <- DESeqDataSetFromTximport(txi = txi1, colData=samp.info, design = ~ 0 + Pop) 
#Remove low reads
dds <- dds[rowMeans(counts(dds))>=min,]
#Remove ribosomes
dds <- dds[!names(dds)%in%ribosomes,] 
################################################################################
# vsdata <- vst(dds, blind=FALSE)
# pdf(paste0("figures/wgcna/PCA_LAT.pdf"))
# plotPCA.custom(vsdata, intgroup="Pop", PCx=1, PCy=2) 
# plotPCA.custom(vsdata, intgroup="Pop", PCx=3, PCy=4)  
# plotPCA.custom(vsdata, intgroup="Temp", PCx=1, PCy=2) 
# plotPCA.custom(vsdata, intgroup="Temp", PCx=3, PCy=4) 
# plotPCA.custom(vsdata, intgroup="Genotype", PCx=1, PCy=2) 
# plotPCA.custom(vsdata, intgroup="Genotype", PCx=3, PCy=4) 
# dev.off()
# ################################################################################
# #Correct genotype effect
write.csv(counts(estimateSizeFactors(dds), normalized=TRUE), "scripts/data/lat_gene_expr_normalized_uncorrected.csv")
write.csv(counts(estimateSizeFactors(dds), normalized=FALSE), "scripts/data/lat_gene_expr_raw_uncorrected.csv")

gene_exp_adj <- sva::ComBat_seq(counts = counts(dds, normalized=FALSE), batch = samp.info$Genotype,
                                group = samp.info$Temp)
dds <- DESeqDataSetFromMatrix(gene_exp_adj, colData = samp.info, design = ~ 0 + Pop)
#counts_raw <- counts(estimateSizeFactors(dds), normalized=FALSE)
# 
# vsdata <- vst(dds, blind=FALSE)
# pdf(paste0("figures/wgcna/PCA_LAT_corrected.pdf"))
# plotPCA.custom(vsdata, intgroup="Pop", PCx=1, PCy=2) 
# plotPCA.custom(vsdata, intgroup="Pop", PCx=3, PCy=4)  
# plotPCA.custom(vsdata, intgroup="Temp", PCx=1, PCy=2) 
# plotPCA.custom(vsdata, intgroup="Temp", PCx=3, PCy=4) 
# plotPCA.custom(vsdata, intgroup="Genotype", PCx=1, PCy=2) 
# plotPCA.custom(vsdata, intgroup="Genotype", PCx=3, PCy=4) 
# dev.off()
################################################################################
expression <-t(counts(estimateSizeFactors(dds), normalized=TRUE))
#As recommended, only keep genes that are strongly expressed :
quant <- unlist(lapply(1:ncol(expression), FUN=function(i)( quantile(expression[,i], prob=c(.75), type=1) )))
expression <- expression[,quant>15]
write.csv(t(expression), "scripts/data/lat_gene_expr_normalized_corrected.csv")
#expression <- t(fpkm(dds))
#expression <- log2(expression+1)
#expression[-Inf] <- 0
################################################################################
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
softPower <- 8
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
write.csv(module_df, "scripts/data/Latitudinal_genes_x_module.csv")
count_lat <- counts(estimateSizeFactors(dds), normalized=TRUE)
################################################################################
#https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(expression, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  tidyr::pivot_longer(-treatment) %>%
  mutate(name = gsub("ME", "", name), name = factor(name, levels = module_order))

library("ggplot2")
g1 <- ggplot(mME, aes(x=treatment, y=name, fill=value)) +
  geom_tile() +  theme_bw() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
g1


mME$Strain <- stringr::str_split(mME$treatment, pattern="-",n=3 , simplify = TRUE)[,1]
mME$Strain <- as.factor(mME$Strain)
data <- mME %>% group_by(Strain, name) %>% summarise(vr=mean(value))

ggplot(data, aes(x=Strain, y=name, fill=vr)) +
  geom_tile() +  theme_bw() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


temp17 <- unlist(lapply(1:nrow(mME), FUN=function(i)( mME[i,1] %in% c("PNZT028-IPO94203-17R2","PNZT028-IPO94203-17R3",
                                                                      "PNZT029-LG2012-17R2","PNZT029-LG2012-17R3","PNZT031-GRI002-17R2","PNZT031-GRI002-17R3",
                                                                      "PNZT032-IPO94218-17R2","PNZT032-IPO94218-17R3","PNZT033-IPO6001-17R2","PNZT033-IPO6001-17R3") )))
temp23 <- unlist(lapply(1:nrow(mME), FUN=function(i)( mME[i,1] %in% c("PNZT028-IPO94203-23R2","PNZT028-IPO94203-23R4",
                                                                      "PNZT029-LG2012-23R2","PNZT029-LG2012-23R4",  
                                                                      "PNZT031-GRI002-23R2","PNZT031-GRI002-23R4", 
                                                                      "PNZT032-IPO94218-23R2","PNZT032-IPO94218-23R4",
                                                                      "PNZT033-IPO6001-23R2","PNZT033-IPO6001-23R4" ) )))

mME[temp17,4] <- 17
mME[temp23,4] <- 23
mME[!(temp23 | temp17),4] <- 20
colnames(mME) <- c(colnames(mME[,1:3]),"Temp")
mME$Temp <- as.factor(mME$Temp)
data <- mME %>% group_by(Temp, name) %>% summarise(vr=mean(value))

g2 <- ggplot(data, aes(x=Temp, y=name, fill=vr)) +
  geom_tile() +  theme_bw() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
g2

moduleTraitCor = cor(MEs0[,1:(ncol(MEs0)-1)], as.numeric(samp.info[,2]), use = "p")

#https://www.biostars.org/p/384281/  "you would focus on the modules from which the regression produces p<0.05"
summary(lm(MEs0$MEtan ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.8390196                                   # 14
summary(lm(MEs0$MEbrown ~ as.numeric(samp.info[,2])))$coefficients[,4] #  3.160495e-06 ( corr = 0.738601386)       # 4
summary(lm(MEs0$MEyellow ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.006270073 ( corr = 0.487624022)        # 6
summary(lm(MEs0$MEcyan ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.00650729 ( corr = 0.485720500)           # 7
summary(lm(MEs0$MEgreen ~ as.numeric(samp.info[,2])))$coefficients[,4] # 1.799347e-07  ( corr = 0.792456066)       # 2
summary(lm(MEs0$MEpurple ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.9537741                                # 17
summary(lm(MEs0$MEmagenta ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.9209046                               # 16
summary(lm(MEs0$MEmidnightblue ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.8671156                          # 15
summary(lm(MEs0$MEgreenyellow ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.9602174                           # 18
summary(lm(MEs0$MEblue ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.3328115                                  # 12
summary(lm(MEs0$MEpink ~ as.numeric(samp.info[,2])))$coefficients[,4] # 1.468278e-09  ( corr = -0.857193045)       # 1
summary(lm(MEs0$MEblack ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.2710060                                 # 11
summary(lm(MEs0$MEsalmon ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.06906483 ( corr = -0.336462921)        # 9
summary(lm(MEs0$MElightcyan ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.01628134 ( corr = -0.435031120)     # 8
summary(lm(MEs0$MEgrey60 ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.2631265                                # 10
summary(lm(MEs0$MEred ~ as.numeric(samp.info[,2])))$coefficients[,4] # 6.918041e-05  ( corr = -0.661349463)        # 5
summary(lm(MEs0$MEturquoise ~ as.numeric(samp.info[,2])))$coefficients[,4] # 2.751480e-06   ( corr = -0.741548769) # 3
summary(lm(MEs0$MEgrey ~ as.numeric(samp.info[,2])))$coefficients[,4] # 0.7497234                                  # 13


pdf(file = "figures/wgcna/Module_corr_categories_LAT_unsigned.pdf", wi = 6, he = 6)
g1
g2
dev.off()
################################################################################
#UNCORRECTED
################################################################################
################################################################################
################################################################################
################################################################################
#Data
quantData <- list.files(path = "All_Data/Latitudinal/", full.names=TRUE, recursive = TRUE, pattern = "t_data.ctab")
quantData <- quantData[c(1:12,19:36)]
quantData <- quantData[order(as.character(quantData), method = c("radix"))]
names(quantData) <- str_split(quantData, "/", simplify = TRUE)[,4] #select column corresponding to sample name
#Design info
samp.info <- read.csv("All_Data/Latitudinal/samples_design.csv", sep = ",")
rownames(samp.info) <- samp.info$X
samp.info <- samp.info[c(1:12,19:36),2:4]
samp.info[,1:ncol(samp.info)] <- lapply(samp.info[,1:ncol(samp.info)] , factor)
################################################################################
#Import count data
tmp <- read_tsv(quantData[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi1 <- tximport(files=quantData, type = "stringtie", tx2gene = tx2gene, readLength=75) 
dds <- DESeqDataSetFromTximport(txi = txi1, colData=samp.info, design = ~ 0 + Pop) 
#Remove low reads
dds <- dds[rowMeans(counts(dds))>=min,]
#Remove ribosomes
dds <- dds[!names(dds)%in%ribosomes,] 
################################################################################
dds <- DESeq(dds)
expression <-t(counts(estimateSizeFactors(dds), normalized=TRUE))
#As recommended, only keep genes that are strongly expressed :
quant <- unlist(lapply(1:ncol(expression), FUN=function(i)( quantile(expression[,i], prob=c(.75), type=1) )))
expression <- expression[,quant>10]
write.csv(t(expression), "scripts/data/lat_gene_expr_normalized_uncorrected.csv")
#expression <- t(fpkm(dds))
#expression <- log2(expression+1)
#expression[-Inf] <- 0
################################################################################
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
write.csv(module_df, "scripts/data/Latitudinal_genes_x_module_UNCORRECTED.csv")