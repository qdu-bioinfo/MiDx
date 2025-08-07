## ----closeg, message = FALSE--------------------------------------------------
library("Biobase")
library("genefilter")
data(sample.ExpressionSet)
igenes <- c(300,333,355,419) ##the interesting genes
closeg <- genefinder(sample.ExpressionSet, igenes, 10, method="euc", scale="none")
names(closeg)

## ----gene-indexing------------------------------------------------------------
closeg$"31539_r_at"
Nms1 <- featureNames(sample.ExpressionSet)[closeg$"31539_r_at"$indices]
Nms1

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

