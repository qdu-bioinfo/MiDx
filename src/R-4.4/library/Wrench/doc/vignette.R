## ----setup, include=F, message=F----------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

require(metagenomeSeq)
require(edgeR)
require(Wrench)
require(DESeq2)


## ----getPackage, eval=FALSE---------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("Wrench")

## ----eval = FALSE-------------------------------------------------------------
#  BiocManager::install(“HCBravoLab/Wrench”)

## ----Load, message=FALSE------------------------------------------------------
library(Wrench)

## ----warning=FALSE------------------------------------------------------------
#extract count and group information for from the mouse microbiome data in the metagenomeSeq package
data(mouseData)
mouseData

counts <- MRcounts( mouseData, norm=FALSE )  #get the counts
counts[1:10,1:2]

group <- pData(mouseData)$diet #get the group/condition vector
head(group)

#Running wrench with defaults
W <- wrench( counts, condition=group  )
compositionalFactors <- W$ccf
normalizationFactors <- W$nf

head( compositionalFactors ) #one factor for each sample
head( normalizationFactors)  #one factor for each sample


## ----warning=FALSE------------------------------------------------------------

# -- If using metagenomeSeq
normalizedObject <- mouseData  #mouseData is already a metagenomeSeq object 
normFactors(normalizedObject) <- normalizationFactors

# -- If using edgeR, we must pass in the compositional factors
edgerobj <- edgeR::DGEList( counts=counts,
                     group = as.matrix(group),
                     norm.factors=compositionalFactors )

# -- If using DESeq/DESeq2
deseq.obj <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                   DataFrame(group),
                                   ~ group )
deseq.obj
sizeFactors(deseq.obj) <- normalizationFactors


## ----warning=FALSE------------------------------------------------------------
time <- as.numeric(as.character(pData(mouseData)$relativeTime))
time.levs <- cut( time, breaks = c(0, 6, 28, 42, 56, 70) )
overall_group <- paste( group, time.levs ) #merge the time information and the group information together
W <- wrench( counts, condition = overall_group )

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

