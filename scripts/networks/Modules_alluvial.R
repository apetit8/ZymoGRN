library(ggplot2)
library(ggalluvial)
library(gridExtra)

moduleColors_ext <- read.csv("scripts/data/Latitudinal_genes_x_module_UNCORRECTED.csv")[,2:3]
moduleColors_ext$Dataset <- "3.External" 
moduleColors_full <- read.csv("scripts/data/FULLDATA2_genes_x_module.csv")[,2:3]
moduleColors_full$Dataset <- "2.Both"
moduleColors_anc <- read.csv("scripts/data/Ancestors_genes_x_module_corrected_exp_data.csv")[,2:3]
moduleColors_anc$Dataset <- "1.Ancestors"


moduleColors_ext <- moduleColors_ext[moduleColors_ext$gene_id %in% moduleColors_full$gene_id,]
moduleColors_anc <- moduleColors_anc[moduleColors_anc$gene_id %in% moduleColors_full$gene_id,]
df.mod <- rbind(moduleColors_anc, moduleColors_full, moduleColors_ext)
df.mod <- df.mod[order(df.mod$colors),]

gg <- ggplot(df.mod, aes(x = Dataset, stratum = colors, fill=colors, alluvium = gene_id, width=0.1), show.legend = FALSE) +
  scale_fill_manual(values = unique(df.mod$colors)) +
  geom_flow(stat = "alluvium", lode.guidance = "zagzig") +
  geom_stratum(width=0.1) + 
  theme_bw() +
  guides(fill="none") 


cairo_pdf("figures/Alluvion_mod.pdf", width=8, height=5)
gg
dev.off()

png("figures/Alluvion_mod.png", width=8000, height=5000)
gg
dev.off()


df.mod1 <- rbind(moduleColors_anc, moduleColors_full)
df.mod1 <- df.mod1[order(df.mod1$colors),]

gg1 <- ggplot(df.mod1, aes(x = Dataset, stratum = colors, fill=colors, alluvium = gene_id, width=0.1), show.legend = FALSE) +
  scale_fill_manual(values = unique(df.mod1$colors)) +
  geom_flow(stat = "alluvium", lode.guidance = "zagzig") +
  geom_stratum() + 
  theme_bw() +
  guides(fill="none") 
  #ggtitle("Modules")+ 
  #scale_x_discrete(breaks = c(0.25,2), labels = c("Ancestors", "Both"))


df.mod2 <- rbind(moduleColors_full, moduleColors_ext)
df.mod2 <- df.mod2[order(df.mod2$colors),]

gg2 <- ggplot(df.mod2, aes(x = Dataset, stratum = colors, fill=colors, alluvium = gene_id, width=0.1), show.legend = FALSE) +
  scale_fill_manual(values = unique(df.mod2$colors)) +
  geom_flow(stat = "alluvium", lode.guidance = "zagzig") +
  geom_stratum() + 
  theme_bw() +
  guides(fill="none") 
  #ggtitle("Modules")+ 
  #scale_x_discrete(breaks = c(0.25,2), labels = c("Both", "External dataset"))

pdf("figures/Alluvion_mod2.pdf", width=40, height=20)
grid.arrange(
  gg1,gg2,
  ncol = 2,
  nrow = 1,
  widths = c(1,1),
  clip = FALSE
)
dev.off()


