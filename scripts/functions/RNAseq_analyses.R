library(readr)
library(DESeq2)
library(tximport)
library(tximportData)
library(stringr)
library(ggplot2)
################################################################################

#Function to subset by string pattern
subset.pattern <- function(data, patterns, var, keep.found = TRUE, useBytes = TRUE){
  data$y <- grepl(pattern = paste0(patterns, collapse="|"), x = data[, var],
                  useBytes = useBytes)
  subdata <- subset(data, y == keep.found)
  subdata$y <- NULL
  subdata
}


plotPCA.custom <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, PCx=1, PCy=2) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, PCx], PC2 = pca$x[, PCy], group = group, intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[PCx:PCy]
    return(d)
  }
  ggplot(data=d, aes(x=PC1, y=PC2, color=group)) +
    geom_point(size=3) + 
    xlab(paste0("PC",PCx,": ", round(percentVar[PCx] * 100), "% variance")) +
    ylab(paste0("PC",PCy,": ", round(percentVar[PCy] * 100), "% variance")) +
    coord_fixed()+theme_bw()
}


