## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
library(BiocParallel)
register(SerialParam())

## ----eval=FALSE---------------------------------------------------------------
# score <- sum(colSums(E[p, ])**2) / length(p)

## ----message=FALSE------------------------------------------------------------
library(GEOquery)
library(limma)

gse200250 <- getGEO("GSE200250", AnnotGPL = TRUE)[[1]]

es <- gse200250
es <- es[, grep("Th2_", es$title)]
es$time <- as.numeric(gsub(" hours", "", es$`time point:ch1`))
es <- es[, order(es$time)]

exprs(es) <- normalizeBetweenArrays(log2(exprs(es)), method="quantile")

es <- es[order(rowMeans(exprs(es)), decreasing=TRUE), ]
es <- es[!duplicated(fData(es)$`Gene ID`), ]
rownames(es) <- fData(es)$`Gene ID`
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

fData(es) <- fData(es)[, c("ID", "Gene ID", "Gene symbol")]

es <- es[head(order(rowMeans(exprs(es)), decreasing=TRUE), 12000), ]
head(exprs(es))

## -----------------------------------------------------------------------------
library(msigdbr)
pathwaysDF <- msigdbr("mouse", collection="H")
pathways <- split(as.character(pathwaysDF$ncbi_gene), pathwaysDF$gs_name)

## ----message=FALSE------------------------------------------------------------
library(fgsea)
set.seed(1)
gesecaRes <- geseca(pathways, exprs(es), minSize = 15, maxSize = 500)

## -----------------------------------------------------------------------------
head(gesecaRes, 10)

## ----fig.width=10, fig.height=4, out.width="100%"-----------------------------
plotCoregulationProfile(pathway=pathways[["HALLMARK_E2F_TARGETS"]], 
                        E=exprs(es), titles = es$title, conditions=es$`time point:ch1`)

## ----fig.width=10, fig.height=4, out.width="100%"-----------------------------
plotCoregulationProfile(pathway=pathways[["HALLMARK_HYPOXIA"]], 
                        E=exprs(es), titles = es$title, conditions=es$`time point:ch1`)



## ----fig.width=10, fig.height=6, out.width="100%"-----------------------------
plotGesecaTable(gesecaRes |> head(10), pathways, E=exprs(es), titles = es$title)

## -----------------------------------------------------------------------------
E <- t(base::scale(t(exprs(es)), scale=FALSE))
pcaRev <- prcomp(E, center=FALSE)
Ered <- pcaRev$x[, 1:10]
dim(Ered)

## -----------------------------------------------------------------------------
set.seed(1)
gesecaResRed <- geseca(pathways, Ered, minSize = 15, maxSize = 500, center=FALSE)
head(gesecaResRed, 10)

## ----fig.width=4, fig.height=4------------------------------------------------
library(ggplot2)
ggplot(data=merge(gesecaRes[, list(pathway, logPvalFull=-log10(pval))],
                  gesecaResRed[, list(pathway, logPvalRed=-log10(pval))])) +
    geom_point(aes(x=logPvalFull, y=logPvalRed)) +
    coord_fixed() + theme_classic()

## -----------------------------------------------------------------------------
library(Seurat)

## ----fig.width=8, fig.height=3.5----------------------------------------------
obj <- readRDS(url("https://alserglab.wustl.edu/files/fgsea/GSE116240.rds"))
obj

newIds <- c("0"="Adventitial MF",
            "3"="Adventitial MF",
            "5"="Adventitial MF",
            "1"="Intimal non-foamy MF",
            "2"="Intimal non-foamy MF",
            "4"="Intimal foamy MF",
            "7"="ISG+ MF",
            "8"="Proliferating cells",
            "9"="T-cells",
            "6"="cDC1",
            "10"="cDC2",
            "11"="Non-immune cells")

obj <- RenameIdents(obj, newIds)

DimPlot(obj) + ggplot2::coord_fixed()

## -----------------------------------------------------------------------------
obj <- SCTransform(obj, verbose = FALSE, variable.features.n = 10000)

## -----------------------------------------------------------------------------
length(VariableFeatures(obj)) # make sure it's a full gene universe of 10000 genes
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,
                rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)

E <- obj@reductions$pca.rev@feature.loadings

## -----------------------------------------------------------------------------
library(msigdbr)

pathwaysDF <- msigdbr("mouse", collection="C2", subcollection = "CP:KEGG_LEGACY")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

# The code above requires msigdbdf package (https://igordot.r-universe.dev/msigdbdf) to be installed to get the full list of pathways.
# For this vignette we use presaved KEGG pathways.
pathways <- readRDS(url("https://alserglab.wustl.edu/files/fgsea/example_kegg_pathways.rds"))

## -----------------------------------------------------------------------------
set.seed(1)
gesecaRes <- geseca(pathways, E, minSize = 5, maxSize = 500, center = FALSE, eps=1e-100)

head(gesecaRes, 10)

## ----fig.width=12, fig.height=7, out.width="100%"-----------------------------
topPathways <- gesecaRes[, pathway] |> head(4)
titles <- sub("KEGG_", "", topPathways)

ps <- plotCoregulationProfileReduction(pathways[topPathways], obj,
                                       title=titles,
                                       reduction="tsne")
cowplot::plot_grid(plotlist=ps[1:4], ncol=2)

## ----fig.width=5, fig.height=3.5, out.width="50%"-----------------------------
plotCoregulationProfileReduction(pathways$KEGG_LYSOSOME, 
                               obj,
                               title=sprintf("KEGG_LYSOSOME (pval=%.2g)",
                                             gesecaRes[match("KEGG_LYSOSOME", pathway), pval]),
                               reduction="tsne")

## ----message=FALSE------------------------------------------------------------
library(Seurat)

obj <- readRDS(url("https://alserglab.wustl.edu/files/fgsea/275_T_seurat.rds"))

## -----------------------------------------------------------------------------
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE, variable.features.n = 10000)

obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,
                rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)

E <- obj@reductions$pca.rev@feature.loadings

## -----------------------------------------------------------------------------
library(msigdbr)
pathwaysDF <- msigdbr("human", collection="H")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

## -----------------------------------------------------------------------------
set.seed(1)
gesecaRes <- geseca(pathways, E, minSize = 15, maxSize = 500, center = FALSE)

head(gesecaRes, 10)

## ----fig.width=10, fig.height=7, out.width="100%"-----------------------------

topPathways <- gesecaRes[, pathway] |> head(4)
titles <- sub("HALLMARK_", "", topPathways)

pt.size.factor <- 1.6


# Starting from Seurat version 5.1.0, the scaling method for spatial plots was changed.
# As mentioned here: https://github.com/satijalab/seurat/issues/9049
# This code provides a workaround to adapt to the new scaling behavior.
if (packageVersion("Seurat") >= "5.1.0"){
    sfactors <- ScaleFactors(obj@images$slice1)
    pt.size.factor <- 1.5 * sfactors$fiducial / sfactors$hires
}

ps <- plotCoregulationProfileSpatial(pathways[topPathways], obj,
                                     title=titles, 
                                     pt.size.factor=pt.size.factor)
cowplot::plot_grid(plotlist=ps, ncol=2)

## ----fig.width=5, fig.height=3.5, out.width="50%"-----------------------------
plotCoregulationProfileSpatial(pathways$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                               obj,
                               pt.size.factor=pt.size.factor,
                               title=sprintf("HALLMARK_OXIDATIVE_PHOSPHORYLATION\n(pval=%.2g)",
                                             gesecaRes[
                                                 match("HALLMARK_OXIDATIVE_PHOSPHORYLATION", pathway),
                                                 pval]))

## ----echo=TRUE----------------------------------------------------------------
sessionInfo()

