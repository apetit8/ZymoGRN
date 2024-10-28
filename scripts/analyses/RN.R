source("scripts/functions/RNAseq_analyses.R")
################################################################################
samp.info <- read.csv("All_Data/samples_design.csv")
samp.info <- samp.info[order(as.character(samp.info$sample_id), method = c("radix")),]
samp.info <- subset(samp.info, ancester==01)
samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
#Counts
raw_counts <- read.csv("scripts/data/counts_01_raw_corrected.csv")
rownames(raw_counts) <- raw_counts$X
raw_counts <- raw_counts[,2:ncol(raw_counts)]
dds_all <- DESeq(DESeqDataSetFromMatrix(raw_counts, colData = samp.info, design = ~ 0 + Pop))
counts_normalized <- counts(estimateSizeFactors(dds_all), normalized=TRUE)
################################################################################
RN.df01 <- data.frame(anc=c(rep(0, nrow(counts_normalized))), Fluct=c(rep(0, nrow(counts_normalized))), S17=c(rep(0, nrow(counts_normalized))), S23=c(rep(0, nrow(counts_normalized))), ExprFluct=c(rep(0, nrow(counts_normalized))), ExprS17=c(rep(0, nrow(counts_normalized))), ExprS23=c(rep(0, nrow(counts_normalized))))


for (i in 1:nrow(counts_normalized)) {
  #lm of ancester at 17°C and 23°C
  RN.df01[i,1] <- lm( c(counts_normalized[i,1], counts_normalized[i,2], counts_normalized[i,41], #17
                        counts_normalized[i,3], counts_normalized[i,4], counts_normalized[i,42], counts_normalized[i,43]) ~ c(rep(17, 3), rep(23, 4))  )$coefficients[2]
  #lm of fluctuating envir at 17°C and 23°C
  RN.df01[i,2] <- lm( c(counts_normalized[i,13], counts_normalized[i,14], counts_normalized[i,37], counts_normalized[i,38], counts_normalized[i,25], counts_normalized[i,26], #17
                        counts_normalized[i,15], counts_normalized[i,16], counts_normalized[i,27], counts_normalized[i,28], counts_normalized[i,39], counts_normalized[i,40]) ~ c(rep(17, 6), rep(23, 6))  )$coefficients[2]
  #lm of S17 at 17°C and 23°C
  RN.df01[i,3] <- lm( c(counts_normalized[i,5], counts_normalized[i,6], counts_normalized[i,17], counts_normalized[i,18], counts_normalized[i,29], counts_normalized[i,30], #17
                        counts_normalized[i,7], counts_normalized[i,8], counts_normalized[i,19], counts_normalized[i,20], counts_normalized[i,31], counts_normalized[i,32]) ~ c(rep(17, 6), rep(23, 6))  )$coefficients[2]
  #lm of S23 at 17°C and 23°C
  RN.df01[i,4] <- lm( c(counts_normalized[i,9], counts_normalized[i,10], counts_normalized[i,21], counts_normalized[i,22], counts_normalized[i,33], counts_normalized[i,34], #17
                        counts_normalized[i,11], counts_normalized[i,12], counts_normalized[i,23], counts_normalized[i,24], counts_normalized[i,35], counts_normalized[i,36]) ~ c(rep(17, 6), rep(23, 6))  )$coefficients[2]
  #Expr of ancester at 17°C and 23°C
  RN.df01[i,5] <- mean( c(counts_normalized[i,33], counts_normalized[i,34], #17
                          counts_normalized[i,31], counts_normalized[i,32], counts_normalized[i,35], counts_normalized[i,36]) )
  #Expr of fluctuating envir at 17°C and 23°C
  RN.df01[i,6] <- mean( c(counts_normalized[i,9], counts_normalized[i,10], counts_normalized[i,20], counts_normalized[i,27], counts_normalized[i,28], #17
                          counts_normalized[i,11], counts_normalized[i,12], counts_normalized[i,21], counts_normalized[i,22], counts_normalized[i,29], counts_normalized[i,30])   )
  #Expr of S17 at 17°C and 23°C
  RN.df01[i,7] <- mean( c(counts_normalized[i,1], counts_normalized[i,2], counts_normalized[i,13], counts_normalized[i,14], #17
                          counts_normalized[i,3], counts_normalized[i,4], counts_normalized[i,15])   )
  #Expr of S23 at 17°C and 23°C
  RN.df01[i,8] <- mean( c(counts_normalized[i,5], counts_normalized[i,6], counts_normalized[i,16], counts_normalized[i,17], counts_normalized[i,29], counts_normalized[i,30], #17
                          counts_normalized[i,7], counts_normalized[i,8], counts_normalized[i,18], counts_normalized[i,19], counts_normalized[i,25], counts_normalized[i,26]) )
  
}



# Get the padj of every DF analyses, keep genes when found at least once with a padj <= 0.05
RN.df01$padjanc <- as.data.frame(results(dds_all, contrast = list(c("PopA01_a_23","Pop01A_b_23_n"), c("PopA01_a_17","Pop01A_b_17_n")), listValues=c(1/2, -1/2)))$padj
RN.df01$padjFluct <- as.data.frame(results(dds_all, contrast = list(c("Pop01F_a_23_n", "Pop01F_b_23", "Pop01F_c_23"),c("Pop01F_a_17_n", "Pop01F_b_17","Pop01F_c_17")), listValues=c(1/3, -1/3)))$padj
RN.df01$padjS17 <- as.data.frame(results(dds_all, contrast = list(c("Pop01S17_a_23_n", "Pop01S17_b_23", "Pop01S17_c_23_n"), c("Pop01S17_a_17_n","Pop01S17_b_17","Pop01S17_c_17_n")), listValues=c(1/3, -1/3)))$padj
RN.df01$padjS23 <- as.data.frame(results(dds_all, contrast = list(c("Pop01S23_a_23", "Pop01S23_b_23_n","Pop01S23_c_23_n"), c("Pop01S23_a_17","Pop01S23_b_17_n","Pop01S23_c_17_n")), listValues=c(1/3, -1/3)))$padj

plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05  )[,1:2])
plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05  )[,c(1,3)])
plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05  )[,c(1,4)])

plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05  )[,1:2])
plot(subset(RN.df01, padjanc <= 0.05 | padjS17 <= 0.05  )[,c(1,3)])
plot(subset(RN.df01, padjanc <= 0.05 | padjS23 <= 0.05  )[,c(1,4)])

boxplot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05  )[,1:2], ylim=c(-1000,1000))
boxplot(subset(RN.df01, padjanc <= 0.05 | padjS17 <= 0.05  )[,c(1,3)], ylim=c(-1000,1000))
boxplot(subset(RN.df01, padjanc <= 0.05 | padjS23 <= 0.05  )[,c(1,4)], ylim=c(-1000,1000))


#All in one boxplot
boxplot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,1:4], ylim=c(-1000,1000))
boxplot(abs(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,1:4]), ylim=c(0,2500))

#All genes
boxplot(RN.df01[,1:4], ylim=c(-1000,1000))
boxplot(abs(RN.df01[,1:4]), ylim=c(0,2500))

#Gene expression / RN
plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,5],abs(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,1]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,6],abs(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,2]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,7],abs(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,3]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,8],abs(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,4]), log = "y", xlim=c(0,30000), ylim=c(1,50000))


plot(RN.df01[,5],abs(RN.df01[,1]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(RN.df01[,6],abs(RN.df01[,2]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(RN.df01[,7],abs(RN.df01[,3]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(RN.df01[,8],abs(RN.df01[,4]), log = "y", xlim=c(0,30000), ylim=c(1,50000))

#Boxplot of expression: nothing special
boxplot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,5:8], ylim=c(0,20000))

################################################################################
#01 by replicate rather than treatment
################################################################################
RN.rep.df01 <- data.frame(anc=c(rep(0, nrow(counts_normalized))),
                          Fluct_a=c(rep(0, nrow(counts_normalized))), Fluct_b=c(rep(0, nrow(counts_normalized))), Fluct_c=c(rep(0, nrow(counts_normalized))),
                          S17a=c(rep(0, nrow(counts_normalized))),S17b=c(rep(0, nrow(counts_normalized))),S17c=c(rep(0, nrow(counts_normalized))),
                          S23a=c(rep(0, nrow(counts_normalized))), S23b=c(rep(0, nrow(counts_normalized))), S23c=c(rep(0, nrow(counts_normalized))))

for (i in 1:nrow(counts_normalized)) {
  #lm of ancester at 17°C and 23°C
  RN.rep.df01[i,1] <- lm( c(counts_normalized[i,1], counts_normalized[i,2], counts_normalized[i,41], #17
                            counts_normalized[i,3], counts_normalized[i,4], counts_normalized[i,42], counts_normalized[i,43]) ~ c(rep(17, 3), rep(23, 4))  )$coefficients[2]
  #lm of fluctuating envir at 17°C and 23°C
  RN.rep.df01[i,2] <- lm( c(counts_normalized[i,13], counts_normalized[i,14], #17
                            counts_normalized[i,15], counts_normalized[i,16]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df01[i,3] <- lm( c( counts_normalized[i,25], counts_normalized[i,26], #17
                             counts_normalized[i,27], counts_normalized[i,28]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df01[i,4] <- lm( c(counts_normalized[i,37], counts_normalized[i,38],  #17
                            counts_normalized[i,39], counts_normalized[i,40]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  #lm of S17 at 17°C and 23°C
  RN.rep.df01[i,5] <- lm( c(counts_normalized[i,5], counts_normalized[i,6], #17
                            counts_normalized[i,7], counts_normalized[i,8]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df01[i,6] <- lm( c(counts_normalized[i,17], counts_normalized[i,18], #17
                            counts_normalized[i,19], counts_normalized[i,20]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df01[i,7] <- lm( c(counts_normalized[i,29], counts_normalized[i,30], #17
                            counts_normalized[i,31], counts_normalized[i,32]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  #lm of S23 at 17°C and 23°C
  RN.rep.df01[i,8] <- lm( c(counts_normalized[i,9], counts_normalized[i,10], #17
                            counts_normalized[i,11], counts_normalized[i,12]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df01[i,9] <- lm( c(counts_normalized[i,21], counts_normalized[i,22], #17
                            counts_normalized[i,23], counts_normalized[i,24]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df01[i,10] <- lm( c(counts_normalized[i,33], counts_normalized[i,34], #17
                             counts_normalized[i,35], counts_normalized[i,36]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
}

# Get the padj of every DF analyses, keep genes when found at least once with a padj <= 0.05
RN.rep.df01$padjanc <- as.data.frame(results(dds_all, contrast = list(c("PopA01_a_23","Pop01A_b_23_n"), c("PopA01_a_17","Pop01A_b_17_n")), listValues=c(1/2, -1/2)))$padj
RN.rep.df01$padjFluct_a <- as.data.frame(results(dds_all, contrast = list(c("Pop01F_a_23_n"),c("Pop01F_a_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjFluct_b <- as.data.frame(results(dds_all, contrast = list(c("Pop01F_b_23"),c("Pop01F_b_17")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjFluct_c <- as.data.frame(results(dds_all, contrast = list(c("Pop01F_c_23"),c("Pop01F_c_17")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjS17a <- as.data.frame(results(dds_all, contrast = list(c("Pop01S17_a_23_n"), c("Pop01S17_a_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjS17b <- as.data.frame(results(dds_all, contrast = list(c("Pop01S17_b_23"), c("Pop01S17_b_17")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjS17c <- as.data.frame(results(dds_all, contrast = list(c("Pop01S17_c_23_n"), c("Pop01S17_c_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjS23a <- as.data.frame(results(dds_all, contrast = list(c("Pop01S23_a_23"), c("Pop01S23_a_17")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjS23b <- as.data.frame(results(dds_all, contrast = list(c("Pop01S23_b_23_n"), c("Pop01S23_b_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df01$padjS23c <- as.data.frame(results(dds_all, contrast = list(c("Pop01S23_c_23_n"), c("Pop01S23_c_17_n")), listValues=c(1/1, -1/1)))$padj

#All in one boxplot
boxplot(subset(RN.rep.df01, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS17c <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:10], ylim=c(-1000,1000))
boxplot(abs(subset(RN.rep.df01, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS17c <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:10]), ylim=c(0,2500))


#For every gene, count how many times plast changed times.
for (i in 1:nrow(RN.rep.df01)) {
  RN.rep.df01[i,21] <- sum(sign(RN.rep.df01[i,1]) != sign(RN.rep.df01[i,2:10]) )
}

hist(RN.rep.df01[,21], seq(-0.5,9.5, 1))
hist(subset(RN.rep.df01, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,21], seq(-0.5,9.5, 1))
hist(subset(RN.rep.df01, abs(anc) <= 5000)$anc, seq(-5000,5000, 100))

#Do same dataset BUT with random RN and see if the distribution is different ?? How to "random" RN ?
#Have to keep same distrib so could pick randomly existing values
RN.rep.df01.rdm <- RN.rep.df01
RN.rep.df01.rdm[,1] <- sample(unlist(RN.rep.df01[,1]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,2] <- sample(unlist(RN.rep.df01[,2]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,3] <- sample(unlist(RN.rep.df01[,3]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,4] <- sample(unlist(RN.rep.df01[,4]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,5] <- sample(unlist(RN.rep.df01[,5]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,6] <- sample(unlist(RN.rep.df01[,6]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,7] <- sample(unlist(RN.rep.df01[,7]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,8] <- sample(unlist(RN.rep.df01[,8]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,9] <- sample(unlist(RN.rep.df01[,9]), size=nrow(RN.rep.df01))
RN.rep.df01.rdm[,10] <- sample(unlist(RN.rep.df01[,10]), size=nrow(RN.rep.df01)) 

#For every gene, count how many times plast changed times.
for (i in 1:nrow(RN.rep.df01.rdm)) {
  RN.rep.df01.rdm[i,21] <- sum(sign(RN.rep.df01.rdm[i,1]) != sign(RN.rep.df01.rdm[i,2:9]) )
}
hist(RN.rep.df01.rdm[,21], seq(-0.5,9.5, 1), main = "Random RN")


pdf("figures/RN_distrib_01.pdf", width = 5, height = 5)
boxplot(subset(RN.df01, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,5:8], ylim=c(0,10000), main="Gene expr distrib", col=c("grey","purple","skyblue","orangered"))
boxplot(RN.df01[,1:4], ylim=c(-500,500), main="RN", col=c("grey","purple","skyblue","orangered"))
boxplot(abs(RN.df01[,1:4]), ylim=c(0,600), main="abs RN", col=c("grey","purple","skyblue","orangered"))
boxplot(subset(RN.rep.df01, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:10], main="RN DEG", ylim=c(-500,500), col=c("grey","purple","purple","purple","skyblue","skyblue","skyblue","orangered","orangered","orangered"))
boxplot(RN.rep.df01[,1:10], main="RN all genes", ylim=c(-500,500), col=c("grey","purple","purple","purple","skyblue","skyblue","skyblue","orangered","orangered","orangered"))
boxplot(abs(subset(RN.rep.df01, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:10]), main="abs RN DEG", ylim=c(0,800), col=c("grey","purple","purple","purple","skyblue","skyblue","skyblue","orangered","orangered","orangered"))
boxplot(abs(RN.rep.df01[,1:10]), main="abs RN all genes", ylim=c(0,800), col=c("grey","purple","purple","purple","skyblue","skyblue","skyblue","orangered","orangered","orangered"))
hist(subset(RN.rep.df01, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,21], seq(-0.5,9.5, 1), xlab = "Number of RN sign changes between Anc & Evolved", ylab="Number of gene", main = "01 RN")
hist(RN.rep.df01.rdm[,21], seq(-0.5,9.5, 1), xlab = "Number of RN sign changes between Anc & Evolved", ylab="Number of gene", main = "Random RN")
dev.off()

################################################################################
#MGGP44
################################################################################
#MGGP44
################################################################################
samp.info <- read.csv("All_Data/samples_design.csv")
samp.info <- samp.info[order(as.character(samp.info$sample_id), method = c("radix")),]
samp.info <- subset(samp.info, ancester==44)
samp.info[,2:ncol(samp.info)] <- lapply(samp.info[,2:ncol(samp.info)] , factor)
#Counts
raw_counts <- read.csv("scripts/data/counts_44_raw_corrected.csv")
rownames(raw_counts) <- raw_counts$X
raw_counts <- raw_counts[,2:ncol(raw_counts)]
dds_all <- DESeq(DESeqDataSetFromMatrix(raw_counts, colData = samp.info, design = ~ 0 + Pop))
counts_normalized <- counts(estimateSizeFactors(dds_all), normalized=TRUE)
################################################################################
RN.df44 <- data.frame(anc=c(rep(0, nrow(counts_normalized))), Fluct=c(rep(0, nrow(counts_normalized))), S17=c(rep(0, nrow(counts_normalized))), S23=c(rep(0, nrow(counts_normalized))),
                      EXPR_anc=c(rep(0, nrow(counts_normalized))), EXPR_Fluct=c(rep(0, nrow(counts_normalized))), EXPR_S17=c(rep(0, nrow(counts_normalized))), EXPR_S23=c(rep(0, nrow(counts_normalized))))


for (i in 1:nrow(counts_normalized)) {
  #lm of ancester at 17°C and 23°C
  RN.df44[i,1] <- lm( c(counts_normalized[i,33], counts_normalized[i,34], #17
                        counts_normalized[i,31], counts_normalized[i,32], counts_normalized[i,35], counts_normalized[i,36]) ~ c(rep(17, 2), rep(23, 4))  )$coefficients[2]
  #lm of Fluct at 17°C and 23°C
  RN.df44[i,2] <- lm( c(counts_normalized[i,9], counts_normalized[i,10], counts_normalized[i,20], counts_normalized[i,27], counts_normalized[i,28], #17
                        counts_normalized[i,11], counts_normalized[i,12], counts_normalized[i,21], counts_normalized[i,22], counts_normalized[i,29], counts_normalized[i,30]) ~ c(rep(17, 5), rep(23, 6))  )$coefficients[2]
  #lm of S17 envir at 17°C and 23°C
  RN.df44[i,3] <- lm( c(counts_normalized[i,1], counts_normalized[i,2], counts_normalized[i,13], counts_normalized[i,14], #17
                        counts_normalized[i,3], counts_normalized[i,4], counts_normalized[i,15]) ~ c(rep(17, 4), rep(23, 3))  )$coefficients[2]
  #lm of S23 at 17°C and 23°C
  RN.df44[i,4] <- lm( c(counts_normalized[i,5], counts_normalized[i,6], counts_normalized[i,16], counts_normalized[i,17], counts_normalized[i,23], counts_normalized[i,24], #17
                        counts_normalized[i,7], counts_normalized[i,8], counts_normalized[i,18], counts_normalized[i,19], counts_normalized[i,25], counts_normalized[i,26]) ~ c(rep(17, 6), rep(23, 6))  )$coefficients[2]
  #Expr of ancester at 17°C and 23°C
  RN.df44[i,5] <- mean( c(counts_normalized[i,33], counts_normalized[i,34], #17
                          counts_normalized[i,31], counts_normalized[i,32], counts_normalized[i,35], counts_normalized[i,36]) )
  #Expr of fluctuating envir at 17°C and 23°C
  RN.df44[i,6] <- mean( c(counts_normalized[i,9], counts_normalized[i,10], counts_normalized[i,20], counts_normalized[i,23], counts_normalized[i,24], #17
                          counts_normalized[i,11], counts_normalized[i,12], counts_normalized[i,21], counts_normalized[i,22], counts_normalized[i,25], counts_normalized[i,26])   )
  #Expr of S17 at 17°C and 23°C
  RN.df44[i,7] <- mean( c(counts_normalized[i,1], counts_normalized[i,2], counts_normalized[i,13], counts_normalized[i,14], #17
                          counts_normalized[i,3], counts_normalized[i,4], counts_normalized[i,15])   )
  #Expr of S23 at 17°C and 23°C
  RN.df44[i,8] <- mean( c(counts_normalized[i,5], counts_normalized[i,6], counts_normalized[i,16], counts_normalized[i,17], counts_normalized[i,29], counts_normalized[i,30], #17
                          counts_normalized[i,7], counts_normalized[i,8], counts_normalized[i,18], counts_normalized[i,19], counts_normalized[i,25], counts_normalized[i,26]) )
}


# Get the padj of every DF analyses, keep genes when found at least once with a padj <= 0.05
RN.df44$padjanc <- as.data.frame(results(dds_all, contrast = list(c("Pop44A_a_23", "Pop44A_b_23_n"), c("Pop44A_a_17")), listValues=c(1/2, -1/1)))$padj
RN.df44$padjFluct <- as.data.frame(results(dds_all, contrast = list(c("Pop44F_a_23", "Pop44F_b_23_n", "Pop44F_c_23"),c("Pop44F_a_17", "Pop44F_b_17_n","Pop44F_c_17")), listValues=c(1/3, -1/3)))$padj
RN.df44$padjS17 <- as.data.frame(results(dds_all, contrast = list(c("Pop44S17_a_23", "Pop44S17_b_23_n"), c("Pop44S17_a_17","Pop44S17_b_17_n")), listValues=c(1/2, -1/2)))$padj
RN.df44$padjS23 <- as.data.frame(results(dds_all, contrast = list(c("Pop44S23_a_23_n", "Pop44S23_b_23_n","Pop44S23_c_23"), c("Pop44S23_a_17_n","Pop44S23_b_17_n","Pop44S23_c_17")), listValues=c(1/3, -1/3)))$padj

plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05  )[,1:2])
plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05  )[,c(1,3)])
plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05  )[,c(1,4)])

plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05  )[,1:2])
plot(subset(RN.df44, padjanc <= 0.05 | padjS17 <= 0.05  )[,c(1,3)])
plot(subset(RN.df44, padjanc <= 0.05 | padjS23 <= 0.05  )[,c(1,4)])

boxplot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05  )[,1:2], ylim=c(-1000,1000))
boxplot(subset(RN.df44, padjanc <= 0.05 | padjS17 <= 0.05  )[,c(1,3)], ylim=c(-1000,1000))
boxplot(subset(RN.df44, padjanc <= 0.05 | padjS23 <= 0.05  )[,c(1,4)], ylim=c(-1000,1000))


#All in one boxplot
boxplot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,1:4], ylim=c(-1000,1000))
boxplot(abs(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,1:4]), ylim=c(0,2500))

#All genes
boxplot(RN.df44[,1:4], ylim=c(-1000,1000))
boxplot(abs(RN.df44[,1:4]), ylim=c(0,2500))

#Gene expression / RN
plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,5],abs(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,1]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,6],abs(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,2]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,7],abs(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,3]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,8],abs(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,4]), log = "y", xlim=c(0,30000), ylim=c(1,50000))


plot(RN.df44[,5],abs(RN.df44[,1]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(RN.df44[,6],abs(RN.df44[,2]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(RN.df44[,7],abs(RN.df44[,3]), log = "y", xlim=c(0,30000), ylim=c(1,50000))
plot(RN.df44[,8],abs(RN.df44[,4]), log = "y", xlim=c(0,30000), ylim=c(1,50000))

#Boxplot of expression: nothing special
boxplot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,5:8], ylim=c(0,5000))

################################################################################
#44 by replicate rather than treatment
################################################################################
RN.rep.df44 <- data.frame(anc=c(rep(0, nrow(counts_normalized))),
                          Fluct_a=c(rep(0, nrow(counts_normalized))), Fluct_b=c(rep(0, nrow(counts_normalized))), Fluct_c=c(rep(0, nrow(counts_normalized))),
                          S17a=c(rep(0, nrow(counts_normalized))),S17b=c(rep(0, nrow(counts_normalized))),
                          S23a=c(rep(0, nrow(counts_normalized))), S23b=c(rep(0, nrow(counts_normalized))), S23c=c(rep(0, nrow(counts_normalized))))

colnames(counts_normalized)
for (i in 1:nrow(counts_normalized)) {
  #lm of ancester at 17°C and 23°C
  RN.rep.df44[i,1] <- lm( c(counts_normalized[i,33], counts_normalized[i,34], #17
                            counts_normalized[i,31], counts_normalized[i,32], counts_normalized[i,35], counts_normalized[i,36]) ~ c(rep(17, 2), rep(23, 4))  )$coefficients[2]
  #lm of Fluct at 17°C and 23°C
  RN.rep.df44[i,2] <- lm( c(counts_normalized[i,9], counts_normalized[i,10], #17
                            counts_normalized[i,11], counts_normalized[i,12]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df44[i,3] <- lm( c(counts_normalized[i,20],  #17
                            counts_normalized[i,21], counts_normalized[i,22]) ~ c(rep(17, 1), rep(23, 2))  )$coefficients[2]
  RN.rep.df44[i,4] <- lm( c( counts_normalized[i,27], counts_normalized[i,28], #17
                             counts_normalized[i,29], counts_normalized[i,30]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  #lm of S17 envir at 17°C and 23°C
  RN.rep.df44[i,5] <- lm( c(counts_normalized[i,1], counts_normalized[i,2], #17
                            counts_normalized[i,3], counts_normalized[i,4]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df44[i,6] <- lm( c( counts_normalized[i,13], counts_normalized[i,14], #17
                             counts_normalized[i,15]) ~ c(rep(1, 2), rep(2, 1))  )$coefficients[2]
  #lm of S23 at 17°C and 23°C
  RN.rep.df44[i,7] <- lm( c(counts_normalized[i,5], counts_normalized[i,6],  #17
                            counts_normalized[i,7], counts_normalized[i,8]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df44[i,8] <- lm( c(counts_normalized[i,16], counts_normalized[i,17], #17
                            counts_normalized[i,18], counts_normalized[i,19]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
  RN.rep.df44[i,9] <- lm( c( counts_normalized[i,23], counts_normalized[i,24], #17
                             counts_normalized[i,25], counts_normalized[i,26]) ~ c(rep(17, 2), rep(23, 2))  )$coefficients[2]
}
# Get the padj of every DF analyses, keep genes when found at least once with a padj <= 0.05
RN.rep.df44$padjanc <- as.data.frame(results(dds_all, contrast = list(c("Pop44A_a_23", "Pop44A_b_23_n"), c("Pop44A_a_17")), listValues=c(1/2, -1/1)))$padj
RN.rep.df44$padjFluct_a <- as.data.frame(results(dds_all, contrast = list(c("Pop44F_a_23"),c("Pop44F_a_17")), listValues=c(1/1, -1/1)))$padj
RN.rep.df44$padjFluct_b <- as.data.frame(results(dds_all, contrast = list(c("Pop44F_b_23_n"),c("Pop44F_b_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df44$padjFluct_c <- as.data.frame(results(dds_all, contrast = list(c("Pop44F_c_23"),c("Pop44F_c_17")), listValues=c(1/1, -1/1)))$padj
RN.rep.df44$padjS17a <- as.data.frame(results(dds_all, contrast = list(c("Pop44S17_a_23"), c("Pop44S17_a_17")), listValues=c(1/1, -1/1)))$padj
RN.rep.df44$padjS17b <- as.data.frame(results(dds_all, contrast = list(c("Pop44S17_b_23_n"), c("Pop44S17_b_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df44$padjS23a <- as.data.frame(results(dds_all, contrast = list(c("Pop44S23_a_23_n"), c("Pop44S23_a_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df44$padjS23b <- as.data.frame(results(dds_all, contrast = list(c("Pop44S23_b_23_n"), c("Pop44S23_b_17_n")), listValues=c(1/1, -1/1)))$padj
RN.rep.df44$padjS23c <- as.data.frame(results(dds_all, contrast = list(c("Pop44S23_c_23"), c("Pop44S23_c_17")), listValues=c(1/1, -1/1)))$padj

#Test
# dds_all <- DESeq(DESeqDataSetFromMatrix(raw_counts, colData = samp.info, design = ~ Lineage * temp))
# resultsNames(dds_all)
# RN.rep.df44$padjS23c <- as.data.frame(results(dds_all, contrast = list(c("Pop_44A_a_23_vs_44A_a_17","Pop_44A_b_23_n_vs_44A_a_17"), c("Pop44S23_c_17")), listValues=c(1/1, -1/1)      ) )$padj

#All in one boxplot
boxplot(subset(RN.rep.df44, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:9], ylim=c(-1000,1000))
boxplot(abs(subset(RN.rep.df44, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:9]), ylim=c(0,2500))

#For every gene, count how many times plast changed times.
for (i in 1:nrow(RN.rep.df44)) {
  RN.rep.df44[i,19] <- sum(sign(RN.rep.df44[i,1]) != sign(RN.rep.df44[i,2:9]) )
}

hist(RN.rep.df44[,19], seq(-0.5,9.5, 1))
hist(subset(RN.rep.df44, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,19], seq(-0.5,9.5, 1))
hist(subset(RN.rep.df44, abs(anc) <= 5000)$anc, seq(-5000,5000, 100))

#Do same dataset BUT with random RN and see if the distribution is different ?? How to "random" RN ?
#Have to keep same distrib so could pick randomly existing values
RN.rep.df44.rdm <- RN.rep.df44
RN.rep.df44.rdm[,1] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,2] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,3] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,4] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,5] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,6] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,7] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,8] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))
RN.rep.df44.rdm[,9] <- sample(unlist(RN.rep.df44[,1:9]), size=nrow(RN.rep.df44))

#For every gene, count how many times plast changed times.
for (i in 1:nrow(RN.rep.df44.rdm)) {
  RN.rep.df44.rdm[i,19] <- sum(sign(RN.rep.df44.rdm[i,1]) != sign(RN.rep.df44.rdm[i,2:9]) )
}
hist(RN.rep.df44.rdm[,19], seq(-0.5,9.5, 1), main = "Random RN")


pdf("figures/RN_distrib_44.pdf", width = 5, height = 5)
boxplot(subset(RN.df44, padjanc <= 0.05 | padjFluct <= 0.05 | padjS17 <= 0.05 | padjS23 <= 0.05)[,5:8], ylim=c(0,10000), main="Gene expr distrib", col=c("grey","purple","skyblue","orangered"))
boxplot(RN.df44[,1:4], ylim=c(-500,500), main="RN", col=c("grey","purple","skyblue","orangered"))
boxplot(abs(RN.df44[,1:4]), ylim=c(0,1000), main="abs RN", col=c("grey","purple","skyblue","orangered"))
boxplot(subset(RN.rep.df44, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:9], main="RN DEG", ylim=c(-500,500), col=c("grey","purple","purple","purple","skyblue","skyblue","orangered","orangered","orangered"))
boxplot(RN.rep.df44[,1:9], main="RN all genes", ylim=c(-500,500), col=c("grey","purple","purple","purple","skyblue","skyblue","orangered","orangered","orangered"))
boxplot(abs(subset(RN.rep.df44, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,1:9]), main="abs RN DEG", ylim=c(0,1000), col=c("grey","purple","purple","purple","skyblue","skyblue","orangered","orangered","orangered"))
boxplot(abs(RN.rep.df44[,1:9]), main="abs RN all genes", ylim=c(0,1000), col=c("grey","purple","purple","purple","skyblue","skyblue","orangered","orangered","orangered"))
hist(subset(RN.rep.df44, padjanc <= 0.05 | padjFluct_a <= 0.05 | padjFluct_b <= 0.05 | padjFluct_c <= 0.05 | padjS17a <= 0.05 | padjS17b <= 0.05 | padjS23a <= 0.05 | padjS23b <= 0.05 | padjS23c <= 0.05)[,19], seq(-0.5,9.5, 1), xlab = "Number of RN sign changes between Anc & Evolved", ylab="Number of gene", main = "44 RN")
hist(RN.rep.df44.rdm[,19], seq(-0.5,9.5, 1), xlab = "Number of RN sign changes between Anc & Evolved", ylab="Number of gene", main = "Random RN")
dev.off()

