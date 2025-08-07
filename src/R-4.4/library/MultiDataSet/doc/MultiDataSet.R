## ----Load libraries, message=FALSE, warning = FALSE---------------------------
library(MultiDataSet)
library(brgedata)
library(GenomicRanges)

## ----New Multi----------------------------------------------------------------
multi <- createMultiDataSet()
multi

## ----Names empty Multi--------------------------------------------------------
names(multi)
length(names(multi))

## ----add_eset 10--------------------------------------------------------------
data("brge_gexp")
brge_gexp
multi2 <- add_eset(multi, brge_gexp, dataset.type = "expression")
multi2
multi

## ----add_eset 2---------------------------------------------------------------
multi2 <- add_eset(multi2, brge_gexp, dataset.type = "expression", dataset.name = "new")
multi2

## ----add_eset add id to eset--------------------------------------------------
brge_gexp2 <- brge_gexp
brge_gexp2$id <- 1:100
multi2 <- add_eset(multi, brge_gexp2, dataset.type = "expression")
multi2

## ----add_eset overwrite, error=TRUE-------------------------------------------
brge_gexp2 <- brge_gexp[, 1:10]
multi2 <- add_eset(multi, brge_gexp, dataset.type = "expression", warnings = FALSE)
multi2
multi2 <- add_eset(multi2, brge_gexp2, dataset.type = "expression", warnings = FALSE, overwrite = FALSE)
multi2
multi2 <- add_eset(multi2, brge_gexp2, dataset.type = "expression", warnings = FALSE, overwrite = TRUE)
multi2

## ----add_eset GRanges---------------------------------------------------------
multi2 <- add_eset(multi, brge_gexp, dataset.type = "expression", warnings = FALSE, GRanges = NA)
multi2

## ----add_rse overwrite--------------------------------------------------------
data("brge_methy")
brge_methy2 <- brge_methy[1:100, ] ### Subset the original set to speed up
multi <- createMultiDataSet()
multi2 <- add_rse(multi, brge_methy, dataset.type = "methylation", warnings = FALSE)
multi2

## ----add genexep--------------------------------------------------------------
multi <- createMultiDataSet()
multi2 <- add_genexp(multi, brge_gexp)
brge_gexp

## ----add genexp---------------------------------------------------------------
multi2 <- add_genexp(multi2, brge_gexp, dataset.name = "2")
multi2

## ----Subsetting intro---------------------------------------------------------
multi <- createMultiDataSet()

# Remove probes without a position before adding the object
multi <- add_methy(multi, brge_methy)
multi <- add_genexp(multi, brge_gexp)
multi <- add_eset(multi, brge_gexp, dataset.type = "test", GRanges = NA)
multi

## ----subset samples-----------------------------------------------------------
samples <- sampleNames(brge_gexp)[76:100]
multi[samples, ]

## ----subset common samples----------------------------------------------------
commonSamples(multi)
length(intersect(sampleNames(brge_gexp), sampleNames(brge_methy)))

## ----subset tables------------------------------------------------------------
multi[, "expression"]
multi[, c("methylation", "test")]

## ----select tables------------------------------------------------------------
multi[["expression"]]
multi[, "expression", drop = TRUE]

## ----Genomic Ranges-----------------------------------------------------------
range <- GRanges("chr17:1-100000")
multi[, , range]

## ----Genomic Ranges 2---------------------------------------------------------
range2 <- GRanges(c("chr17:1-100000", "chr17:1000000-2000000"))
multi[, , range2]

## ----combined-----------------------------------------------------------------
multi[samples, "expression", range]
multi[samples, "methylation", range, drop = TRUE]

## ----Advanced genes-----------------------------------------------------------
subset(multi, genes == "SLC35E2")

## ----Advanced genes 2---------------------------------------------------------
subset(multi, genes %in% c("SLC35E2", "IPO13", "TRPV1"))
subset(multi, genes == "EEF1A1" | genes == "LPP")

## ----Advanced pheno-----------------------------------------------------------
subset(multi, , sex == "Female")

## ----Combined advanced--------------------------------------------------------
subset(multi, genes == "SLC35E2", sex == "Female")

