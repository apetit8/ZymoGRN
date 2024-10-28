library(matrixStats)
library(RColorBrewer)
library(scales)
#####Local modules##############################################################
##### Mean RN difference ; Internal modules
RN.rep.df01 <- read.csv("RN_01.csv")[,c(2:12)]
moduleColors <- read.csv("scripts/data/FULLDATA2_genes_x_module.csv")[,2:3]
df.RN_mod <- merge(RN.rep.df01, moduleColors, by="gene_id")[,2:12]
#
diffmat_01 <- lapply(unique(moduleColors$colors), function(mod){
  df <- subset(df.RN_mod, colors==mod)
  diff <- data.frame(S01Fluct_a=mean(abs(df$S01anc - df$S01Fluct_a)),
                     S01Fluct_b=mean(abs(df$S01anc - df$S01Fluct_b)),
                     S01Fluct_c=mean(abs(df$S01anc - df$S01Fluct_c)),
                     S01S17a=mean(abs(df$S01anc - df$S01S17a)),
                     S01S17b=mean(abs(df$S01anc - df$S01S17b)),
                     S01S17c=mean(abs(df$S01anc - df$S01S17c)),
                     S01S23a=mean(abs(df$S01anc - df$S01S23a)),
                     S01S23b=mean(abs(df$S01anc - df$S01S23b)),
                     S01S23c=mean(abs(df$S01anc - df$S01S23c))
  )
  return(diff)
})


unlist(diffmat_01)
diffmat_mat_01 <- t(matrix(data=unlist(diffmat_01), nrow=9, ncol = length(unique(moduleColors$colors))))
colnames(diffmat_mat_01) <- colnames(df.RN_mod)[2:10]
rownames(diffmat_mat_01) <- unique(moduleColors$colors)


m <- diffmat_mat_01
m1 <- m
m1[m1 <= 0.02] <- 0.02
m1[m1 >= 0.2] <- 0.2

image(1:ncol(m1), 1:nrow(m1), t(m1), col = heat.colors(256), axes = FALSE, xlab = "", ylab="", main="MGGP01")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), rownames(m),las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],3))

####
#####MGGP44
RN.rep.df44 <- read.csv("RN_44.csv")[,c(2:3,5:6,8:11)]
df.RN_mod <- merge(RN.rep.df44, moduleColors, by="gene_id")[,2:9]
#
diffmat_44 <- lapply(unique(moduleColors$colors), function(mod){
  df <- subset(df.RN_mod, colors==mod)
  diff <- data.frame(S44Fluct_a=mean(abs(df$S44anc - df$S44Fluct_a)),
                     S44Fluct_c=mean(abs(df$S44anc - df$S44Fluct_c)),
                     S44S17a=mean(abs(df$S44anc - df$S44S17a)),
                     S44S23a=mean(abs(df$S44anc - df$S44S23a)),
                     S44S23b=mean(abs(df$S44anc - df$S44S23b)),
                     S44S23c=mean(abs(df$S44anc - df$S44S23c))
  )
  return(diff)
})


unlist(diffmat_44)
diffmat_mat_44 <- t(matrix(data=unlist(diffmat_44), nrow=6, ncol = length(unique(moduleColors$colors))))
colnames(diffmat_mat_44) <- colnames(df.RN_mod)[2:7]
rownames(diffmat_mat_44) <- unique(moduleColors$colors)


m <- diffmat_mat_44
m1 <- m
m1[m1 <= 0.02] <- 0.02
m1[m1 >= 0.13] <- 0.13

image(1:ncol(m1), 1:nrow(m1), t(m1), col = heat.colors(256), axes = FALSE, xlab = "", ylab="", main="MGGP44")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), rownames(m),las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],3))


#################
pdf("figures/Heatmaps_RN_diff.pdf", width = 12, height = 6)
par(mar = c(3.7,4, 2,0.5), mfrow=c(1,2), fig=c(0.025,0.6,0.1,1))
###
m <- diffmat_mat_01
m <- cbind(m,rowSums2(m))
colnames(m)[10] <- "SUM"
m1 <- m
m1[m1 <= 0.3] <- 0.3
m1[m1 >= 0.8] <- 0.8
image(1:ncol(m1), 1:nrow(m1), t(m1), col = heat.colors(256), axes = FALSE, xlab = "", ylab="", main="MGGP01")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), rownames(m),las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],2))
###
par(mar = c(3.7,0.5, 2,2), fig=c(0.6,1,0.1,1), new=TRUE)
m <- diffmat_mat_44
m <- cbind(m,rowSums2(m))
colnames(m)[7] <- "SUM"
m1 <- m
m1[m1 <= 0.3] <- 0.3
m1[m1 >= 0.8] <- 0.8

image(1:ncol(m1), 1:nrow(m1), t(m1), col = heat.colors(256), axes = FALSE, xlab = "", ylab="", main="MGGP44")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), labels = FALSE,las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],2))
dev.off()

################################################################################
################################################################################
#Plasticity Intensity
################################################################################
################################################################################
##### Mean RN difference ; External modules
RN.rep.df01 <- read.csv("RN_01.csv")[,c(2:12)]
df.RN_mod <- merge(RN.rep.df01, moduleColors, by="gene_id")[,2:12]
#
RNmatEXT_01 <- lapply(unique(moduleColors$colors), function(mod){
  df <- subset(df.RN_mod, colors==mod)
  diff <- data.frame(Ancestor=mean(abs(df$S01anc)),
                     S01Fluct_a=mean(abs(df$S01Fluct_a)),
                     S01Fluct_b=mean(abs(df$S01Fluct_b)),
                     S01Fluct_c=mean(abs(df$S01Fluct_c)),
                     S01S17a=mean(abs(df$S01S17a)),
                     S01S17b=mean(abs(df$S01S17b)),
                     S01S17c=mean(abs(df$S01S17c)),
                     S01S23a=mean(abs(df$S01S23a)),
                     S01S23b=mean(abs(df$S01S23b)),
                     S01S23c=mean(abs(df$S01S23c))
  )
  return(diff)
})
#names(RNmatEXT_01) <- unique(moduleColors$colors)
#unlist(RNmatEXT_01)
RNmatEXT_mat_01 <- matrix(data=unlist(RNmatEXT_01), byrow = TRUE, ncol = 10 )
colnames(RNmatEXT_mat_01) <- colnames(RNmatEXT_01[[1]])
rownames(RNmatEXT_mat_01) <- unique(moduleColors$colors)


m <- RNmatEXT_mat_01
m1 <- m
m1[m1 <= 0.05] <- 0.05
m1[m1 >= 0.15] <- 0.15

image(1:ncol(m1), 1:nrow(m1), t(m1), col = terrain.colors(256), axes = FALSE, xlab = "", ylab="", main="MGGP01")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), rownames(m),las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],3))

####
#####Local modules##############################################################
##### Mean RN difference ; Experimenta modules
RN.rep.df44 <- read.csv("RN_44.csv")[,c(2:3,5:6,8:11)]
df.RN_mod <- merge(RN.rep.df44, moduleColors, by="gene_id")[,2:9]
#
RNmatEXT_44 <- lapply(unique(moduleColors$colors), function(mod){
  df <- subset(df.RN_mod, colors==mod)
  diff <- data.frame(Ancestor=mean(abs(df$S44anc)),
                     S44Fluct_a=mean(abs(df$S44Fluct_a)),
                     S44Fluct_c=mean(abs(df$S44Fluct_c)),
                     S44S17a=mean(abs(df$S44S17a)),
                     S44S23a=mean(abs(df$S44S23a)),
                     S44S23b=mean(abs(df$S44S23b)),
                     S44S23c=mean(abs(df$S44S23c))
  )
  return(diff)
})

RNmatEXT_mat_44 <- matrix(data=unlist(RNmatEXT_44), byrow = TRUE, ncol = 7 )
colnames(RNmatEXT_mat_44) <- colnames(RNmatEXT_44[[1]])
rownames(RNmatEXT_mat_44) <- unique(moduleColors$colors)


m <- RNmatEXT_mat_44
m1 <- m
m1[m1 <= 0.05] <- 0.05
m1[m1 >= 0.12] <- 0.12

image(1:ncol(m1), 1:nrow(m1), t(m1), col = colorRampPalette(brewer.pal(8, "Blues"))(200), axes = FALSE, xlab = "", ylab="", main="MGGP44")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), rownames(m),las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],3))


#################
pdf("figures/Heatmaps_RN_mean.pdf", width = 12, height = 6)
###
par(mar = c(3.7,4, 2,0.5), mfrow=c(1,2), fig=c(0.025,0.6,0.1,1))
m <- RNmatEXT_mat_01
# m <- cbind(m,rowMeans2(m))
# colnames(m)[11] <- "MEAN"
m1 <- m
m1[m1 <= 0.3] <- 0.3
m1[m1 >= 0.8] <- 0.8
image(1:ncol(m1), 1:nrow(m1), t(m1), col = colorRampPalette(brewer.pal(8, "Blues"))(200), axes = FALSE, xlab = "", ylab="", main="MGGP01")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), rownames(m),las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],2))
###
par(mar = c(3.7,0.5, 2,2), fig=c(0.6,1,0.1,1), new=TRUE)
m <- RNmatEXT_mat_44
# m <- cbind(m,rowMeans2(m))
# colnames(m)[8] <- "MEAN"
m1 <- m
m1[m1 <= 0.3] <- 0.3
m1[m1 >= 0.8] <- 0.8
image(1:ncol(m1), 1:nrow(m1), t(m1), col = colorRampPalette(brewer.pal(8, "Blues"))(200), axes = FALSE, xlab = "", ylab="", main="MGGP44")
axis(1, 1:ncol(m), colnames(m),las=2)
axis(2, 1:nrow(m), labels = FALSE, las=2)
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, round(m[y,x],2))
dev.off()

################################################################################
################################################################################
pdf("figures/Reg_anc_plast.pdf", width = 5, height = 4)
par(mar = c(4,4, 2,0.5), mfrow=c(1,1))

plot(RNmatEXT_mat_01[,1], rowMeans2(diffmat_mat_01), pch=19, col=rownames(RNmatEXT_mat_01), xlab="Mean ancestral RN", ylab="Mean RN difference")
points(RNmatEXT_mat_44[,1], rowMeans2(diffmat_mat_44), pch=17, col=rownames(RNmatEXT_mat_01))
abline(lm(c(rowMeans2(diffmat_mat_01),rowMeans2(diffmat_mat_44))~c(RNmatEXT_mat_01[,1],RNmatEXT_mat_44[,1])  ))
cortt <- cor.test(c(RNmatEXT_mat_01[,1],RNmatEXT_mat_44[,1]),c(rowMeans2(diffmat_mat_01),rowMeans2(diffmat_mat_44)))
text(0.7,0.8, paste0("Corr = ", round(cortt$estimate,2), " ***" ))
legend(0.95, 0.5, legend=c("MGGP01", "MGGP44"),
       pch=c(19, 17), cex=0.8)

dev.off()
################################################################################
################################################################################
RN.rep.df01 <- read.csv("RN_01.csv")[,c(2:12)]
moduleColors <- read.csv("scripts/data/FULLDATA2_genes_x_module.csv")[,2:3]
df.RN_mod01 <- merge(RN.rep.df01, moduleColors, by="gene_id")[,2:12]

#
df01  <- data.frame(S01Fluct_a=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01Fluct_a)),
                 S01Fluct_b=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01Fluct_b)),
                 S01Fluct_c=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01Fluct_c)),
                 S01S17a=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01S17a)),
                 S01S17b=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01S17b)),
                 S01S17c=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01S17c)),
                 S01S23a=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01S23a)),
                 S01S23b=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01S23b)),
                 S01S23c=(abs(df.RN_mod01$S01anc - df.RN_mod01$S01S23c)))

#####MGGP44
RN.rep.df44 <- read.csv("RN_44.csv")[,c(2:3,5:6,8:11)]
df.RN_mod44 <- merge(RN.rep.df44, moduleColors, by="gene_id")[,2:9]
#
df44 <- data.frame(S44Fluct_a=(abs(df.RN_mod44$S44anc - df.RN_mod44$S44Fluct_a)),
                     S44Fluct_c=(abs(df.RN_mod44$S44anc - df.RN_mod44$S44Fluct_c)),
                     S44S17a=(abs(df.RN_mod44$S44anc - df.RN_mod44$S44S17a)),
                     S44S23a=(abs(df.RN_mod44$S44anc - df.RN_mod44$S44S23a)),
                     S44S23b=(abs(df.RN_mod44$S44anc - df.RN_mod44$S44S23b)),
                     S44S23c=(abs(df.RN_mod44$S44anc - df.RN_mod44$S44S23c)))


pdf("figures/Plast_evol_reg.pdf", width = 10, height = 4)
par(mar = c(4,4, 2,0.5), mfrow=c(1,2))

plot(abs(df.RN_mod01$S01anc), rowMeans(df01), xlab="genes ancestral RN", ylab="genes mean RN evolution",pch=19, col=alpha("black", 0.2), main="MGGP01")
abline(lm(rowMeans(df01)~df.RN_mod01$S01anc  ), col="red")
cortt <- cor.test(df.RN_mod01$S01anc,rowMeans(df01))
text(1.1,4.5, paste0("Corr = ", round(cortt$estimate,2), " ***" ))
#
plot(abs(df.RN_mod44$S44anc), rowMeans(df44), xlab="genes ancestral RN", ylab="genes mean RN evolution",pch=19, col=alpha("black", 0.2), main="MGGP44")
abline(lm(rowMeans(df44)~df.RN_mod44$S44anc  ), col="red")
cortt <- cor.test(df.RN_mod44$S44anc,rowMeans(df44))
text(1.1,3.5, paste0("Corr = ", round(cortt$estimate,2), " ***" ))
dev.off()



