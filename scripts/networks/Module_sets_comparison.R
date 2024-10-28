exp_counts <- read.csv("All_Data/genes_count_corrected.csv", sep=",")[,c(1:8, 45:50)]
lat_counts2 <- read.csv("scripts/data/lat_gene_expr_normalized_uncorrected.csv")
#
lat_counts2 <- lat_counts2[lat_counts2$X %in% exp_counts$X,]
exp_counts <- exp_counts[exp_counts$X %in% lat_counts2$X,]
#

################################################################################
################################################################################
################################################################################
moduleColors <- read.csv("scripts/data/Ancestors_genes_x_module_corrected_exp_data.csv")[,2:3]
#
################################################################################
#Lat boxplot
################################################################################

lat_corr2 <- lapply(unique(moduleColors$colors), function(mod){
  df <- log2(lat_counts2[(lat_counts2$X %in% moduleColors[moduleColors[,2]==mod,1]),2:31]+1)
  #df[df==-Inf] <- 0
  return( abs(cor(t(df)) )  )
})
names(lat_corr2) <- unique(moduleColors$colors)


################################################################################
#Exp boxplot
################################################################################
exp_corr <- lapply(unique(moduleColors$colors), function(mod){
  df <- exp_counts[(exp_counts$X %in% moduleColors[moduleColors[,2]==mod,1]),2:13]
  return( abs(cor(t(df)) )  )
})
names(exp_corr) <- unique(moduleColors$colors)

lat_exp_mod_corr2 <- lapply(1:length(unique(moduleColors$colors)), function(i){
  return( c(mean(lat_corr2[[i]]), mean(exp_corr[[i]]))  )
})
cc2 <- matrix(unlist(lat_exp_mod_corr2), ncol=2, byrow = TRUE)
cc2 <- as.data.frame.matrix(cbind(cc2, unique(moduleColors$colors)))

########################################################################
########################################################################
moduleColors <- read.csv("scripts/data/Latitudinal_genes_x_module_UNCORRECTED.csv")[,2:3]
#
################################################################################
#Lat boxplot
################################################################################

lat_corr2 <- lapply(unique(moduleColors$colors), function(mod){
  df <- log2(lat_counts2[(lat_counts2$X %in% moduleColors[moduleColors[,2]==mod,1]),2:31]+1)
  #df[df==-Inf] <- 0
  return( abs(cor(t(df)) )  )
})
names(lat_corr2) <- unique(moduleColors$colors)


################################################################################
#Exp boxplot
################################################################################
exp_corr <- lapply(unique(moduleColors$colors), function(mod){
  df <- exp_counts[(exp_counts$X %in% moduleColors[moduleColors[,2]==mod,1]),2:13]
  return( abs(cor(t(df)) )  )
})
names(exp_corr) <- unique(moduleColors$colors)


lat_exp_mod_corr3 <- lapply(1:length(unique(moduleColors$colors)), function(i){
  return( c(mean(lat_corr2[[i]]), mean(exp_corr[[i]]))  )
})

cc3 <- matrix(unlist(lat_exp_mod_corr3), ncol=2, byrow = TRUE)
cc3 <- as.data.frame.matrix(cbind(cc3, unique(moduleColors$colors)))



########################################################################
########################################################################
########################################################################
moduleColors <- read.csv("scripts/data/FULLDATA2_genes_x_module.csv")[,2:3]
#
################################################################################
#Lat boxplot
################################################################################

lat_corr2 <- lapply(unique(moduleColors$colors), function(mod){
  df <- log2(lat_counts2[(lat_counts2$X %in% moduleColors[moduleColors[,2]==mod,1]),2:31]+1)
  #df[df==-Inf] <- 0
  return( abs(cor(t(df)) )  )
})
names(lat_corr2) <- unique(moduleColors$colors)


################################################################################
#Exp boxplot
################################################################################
exp_corr <- lapply(unique(moduleColors$colors), function(mod){
  df <- exp_counts[(exp_counts$X %in% moduleColors[moduleColors[,2]==mod,1]),2:13]
  return( abs(cor(t(df)) )  )
})
names(exp_corr) <- unique(moduleColors$colors)

lat_exp_mod_corr4 <- lapply(1:length(unique(moduleColors$colors)), function(i){
  return( c(mean(lat_corr2[[i]]), mean(exp_corr[[i]]))  )
})

cc4 <- matrix(unlist(lat_exp_mod_corr4), ncol=2, byrow = TRUE)
cc4 <- as.data.frame.matrix(cbind(cc4, unique(moduleColors$colors)))

########################################################################

pdf("figures/Mod_corr_lm_mean.pdf", width = 4.6, height = 4.6)
par(mar = c(4,4, 2,0.1))

plot(cc2[,1], cc2[,2], ylab="mean abs. correlation in ancestors dataset",
     xlab="mean abs. correlation in external dataset",
     col=cc2[,3], pch=19, main="Modules from Ancestors", xlim=c(0.15,0.52), ylim=c(0.3,0.68)) #
points(cc2[,1], cc2[,2], col="black", cex=1.3)
#abline(lm(as.numeric(cc2[,2]) ~ as.numeric(cc2[,1]) ), col="black")
#text(0.05,0.2,paste0("R2 = ",round((cor(as.numeric(cc2[,1]), as.numeric(cc2[,2]))^2),2)))

#
plot(cc3[,1], cc3[,2], ylab="mean abs. correlation in ancestors dataset",
     xlab="mean abs. correlation in external dataset",
     col=cc3[,3], pch=19, main="Modules from External dataset", xlim=c(0.15,0.52), ylim=c(0.3,0.68))
points(cc3[,1], cc3[,2], col="black", cex=1.3)
#abline(lm(as.numeric(cc3[,2]) ~ as.numeric(cc3[,1]) ), col="black")
#text(0.05,0.12,paste0("R2 = ",round((cor(as.numeric(cc3[,1]), as.numeric(cc3[,2]))^2),2)))
#
plot(cc4[,1], cc4[,2], ylab="mean abs. correlation in ancestors dataset",
     xlab="mean abs. correlation in external dataset",
     col=cc4[,3], pch=19, main="Modules from both dataset", xlim=c(0.15,0.52), ylim=c(0.3,0.68)) #, xlim=c(0,0.35), ylim=c(0,0.25)
points(cc4[,1], cc4[,2], col="black", cex=1.3)
#abline(lm(as.numeric(cc4[,2]) ~ as.numeric(cc4[,1]) ), col="black")
#text(0.05,0.2,paste0("R2 = ",round((cor(as.numeric(cc4[,1]), as.numeric(cc4[,2]))^2),2)))
#
dev.off()

mean(as.numeric(cc4[,2]))
mean(as.numeric(cc3[,2]))
mean(as.numeric(cc2[,2]))
