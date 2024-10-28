################################################################################
moduleColors <- read.csv("scripts/data/FULLDATA2_genes_x_module.csv")[,2:3]
annot <- read.csv("All_Data/Annotations-metadata.csv", sep="")
annot$GOs <- stringr::str_replace_all(annot$GOs, "C:", "; ")
annot$GOs <- stringr::str_replace_all(annot$GOs, "F:", "; ")
annot$GOs <- stringr::str_replace_all(annot$GOs, "P:", "; ")
annot[annot$GOs %in% c(1:100),8] <- "no GO terms"

annot_all <- subset(annot, transcript %in% moduleColors$gene_id   )
tab_annot_all <- table(unlist((stringr::str_split(annot_all$GOs,";", simplify = TRUE))))
#(rank(round(table(annot_all),3))[1:10])/nrow(annot_all)*100
round(tab_annot_all[order(tab_annot_all,decreasing = TRUE)],3)[1:20] #/nrow(annot_all)*100)



#
GO_mod <- lapply(unique(moduleColors$colors), function(mod){
  annot_mod <- subset(annot, transcript %in% subset(moduleColors, colors==mod)$gene_id   )
  tab_annot_mod <- table(unlist((stringr::str_split(annot_mod$GOs,";", simplify = TRUE))))
  tab_annot_mod <- (round(tab_annot_mod[order(tab_annot_mod,decreasing = TRUE)], 3)) #/nrow(annot_mod)*100
  return(c(nrow(annot_mod), sum(grepl("ranscript", annot_mod$Seq..Description, fixed = TRUE)) ,tab_annot_mod))
})
names(GO_mod) <- unique(moduleColors$colors)

#where x is the number of DE genes in the gene set, m is the total number of genes
#in the gene set, n is the total number of genes not in the gene set, and k is the total number of DE genes.
#phyper(27, 412, 10043-412, GO_mod$brown[1], lower.tail = FALSE) #"brown"

GO_ORA <- lapply( 1:length(GO_mod), function(mod){
    GO_ORA_mod <- lapply( names(GO_mod[[mod]])[7: sum((GO_mod[[mod]])> (length(GO_mod[[mod]])*0.02) ) ], function(GO){
      phyper(GO_mod[[mod]][names(GO_mod[[mod]])==GO],
             tab_annot_all[names(tab_annot_all)==GO],
             10043-tab_annot_all[names(tab_annot_all)==GO],
             GO_mod[[mod]][1], lower.tail = FALSE) #"red" GO:0003824
    })
    return(unlist(GO_ORA_mod)[order(unlist(GO_ORA_mod))]) #return(p.adjust(unlist(GO_ORA_mod), method = "bonferroni")[order(p.adjust(unlist(GO_ORA_mod), method = "bonferroni"))])
})
names(GO_ORA) <- names(GO_mod)

tt <- p.adjust(unlist(GO_ORA), method = "fdr")
tt[tt<=0.05]
