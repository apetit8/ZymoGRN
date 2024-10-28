################################################################################
RN.rep.df01 <- read.csv("RN_01.csv")
RN.rep.df44 <- read.csv("RN_44.csv")
moduleColors <- read.csv("scripts/data/FULLDATA2_genes_x_module.csv")[,2:3]
###

df.RN_all <- cbind(RN.rep.df01[,c(2:11)], RN.rep.df44[,c(2:3,5:6,8:(ncol(RN.rep.df44)-1))])
df.RN_all$gene_id <- RN.rep.df01$gene_id
#df.RN_mod <- merge(df.RN_all, moduleColors, by="gene_id")[,2:(ncol(df.RN_all)+1)]
# rownames(df.RN_mod) <-  merge(df.RN_all, moduleColors, by="gene_id")[,1]

RN_black <- subset(df.RN_all, gene_id %in% subset(moduleColors, colors=="black")[,1]  )[,1:(ncol(df.RN_all)-1)]
RN_magenta <- subset(df.RN_all, gene_id %in% subset(moduleColors, colors=="magenta")[,1]  )[,1:(ncol(df.RN_all)-1)]

pdf(paste0("figures/Mod_RN_chosen_modules.pdf"), width = 5, height = 3.5)
par(mar = c(5.2,4,1.2,1))
boxplot(RN_black, main=paste0("Black module, ",nrow(RN_black)," genes"), ylab="Gene reaction norms",xaxt="n", cex=0.8,
        col= c("grey20","pink","pink2","pink3","royalblue1","royalblue","blue3","orangered","red","red3",
               "grey20","pink","pink3","royalblue1","orangered","red","red3"), outline=FALSE, las=2, at = c(0,1.2,2,2.8,4,4.8,5.6,6.8,7.6,8.4,9.6,10.8,11.6,12.8,14,14.8,15.6))
abline(v=9, col="black",lty=3)
axis(1, at = c(0,1.2,2,2.8,4,4.8,5.6,6.8,7.6,8.4,9.6,10.8,11.6,12.8,14,14.8,15.6), las=2, y = par("usr")[3] - 0.45,
     labels = c("Ancestor ", "a", "b", "c", "a", "b", "c", "a", "b", "c", "Ancestor ", "a", "c", "a", "a", "b", "c"))
text(2, 0.5, substitute(paste("MGGP01")), xpd=TRUE, font=3, cex=1)
text(12.5, 0.5, substitute(paste("MGGP44")), xpd=TRUE, font=3, cex=1)
text(2, -0.8, substitute(paste("Fluct.")), xpd=TRUE, font=3, cex=1)
text(4.8, -0.8, substitute(paste("17°C")), xpd=TRUE, font=3, cex=1)
text(7.6, -0.8, substitute(paste("23°C")), xpd=TRUE, font=3, cex=1)
text(11.2, -0.8, substitute(paste("Fluct.")), xpd=TRUE, font=3, cex=1)
text(12.9, -0.8, substitute(paste("17°C")), xpd=TRUE, font=3, cex=1)
text(15, -0.8, substitute(paste("23°C")), xpd=TRUE, font=3, cex=1)


boxplot(RN_magenta, main=paste0("Magenta module, ",nrow(RN_magenta)," genes"), ylab="Gene reaction norms",xaxt="n",cex=0.8,
        col= c("grey20","pink","pink2","pink3","royalblue1","royalblue","blue3","orangered","red","red3",
               "grey20","pink","pink3","royalblue1","orangered","red","red3"), outline=FALSE, las=2, at = c(0,1.2,2,2.8,4,4.8,5.6,6.8,7.6,8.4,9.6,10.8,11.6,12.8,14,14.8,15.6))
abline(v=9, col="black",lty=3)
axis(1, at = c(0,1.2,2,2.8,4,4.8,5.6,6.8,7.6,8.4,9.6,10.8,11.6,12.8,14,14.8,15.6), las=2, y = par("usr")[3] - 0.45,
     labels = c("Ancestor ", "a", "b", "c", "a", "b", "c", "a", "b", "c", "Ancestor ", "a", "c", "a", "a", "b", "c"))
text(2, 0.55, substitute(paste("MGGP01")), xpd=TRUE, font=3, cex=1)
text(12.5, 0.55, substitute(paste("MGGP44")), xpd=TRUE, font=3, cex=1)
text(2, -0.8, substitute(paste("Fluct.")), xpd=TRUE, font=3, cex=1)
text(4.8, -0.8, substitute(paste("17°C")), xpd=TRUE, font=3, cex=1)
text(7.6, -0.8, substitute(paste("23°C")), xpd=TRUE, font=3, cex=1)
text(11.2, -0.8, substitute(paste("Fluct.")), xpd=TRUE, font=3, cex=1)
text(12.9, -0.8, substitute(paste("17°C")), xpd=TRUE, font=3, cex=1)
text(15, -0.8, substitute(paste("23°C")), xpd=TRUE, font=3, cex=1)

dev.off()
