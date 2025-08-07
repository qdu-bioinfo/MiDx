## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width = 6,
                      fig.height = 6,
                      fig.path = 'figures/')

## ----load, message = FALSE----------------------------------------------------
library(ropls)

## ----sacurine, message = FALSE------------------------------------------------
data(sacurine)
names(sacurine)

## ----strF---------------------------------------------------------------------
view(sacurine$dataMatrix)
view(sacurine$sampleMetadata)
view(sacurine$variableMetadata)

## ----pca_code, eval = FALSE---------------------------------------------------
#  sacurine.pca <- opls(sacurine$dataMatrix)

## ----pca_result, echo = FALSE-------------------------------------------------
sacurine.pca <- opls(sacurine$dataMatrix, fig.pdfC = "none")

## ----pca_figure, echo = FALSE, fig.show = 'hold'------------------------------
plot(sacurine.pca)

## ----pca-col------------------------------------------------------------------
genderFc <- sacurine$sampleMetadata[, "gender"]
plot(sacurine.pca,
     typeVc = "x-score",
     parAsColFcVn = genderFc)

## ----pca-col-personalized-----------------------------------------------------
plot(sacurine.pca,
     typeVc = "x-score",
     parAsColFcVn = genderFc,
     parLabVc = as.character(sacurine$sampleMetadata[, "age"]),
     parPaletteVc = c("green4", "magenta"))

## ----plsda--------------------------------------------------------------------
sacurine.plsda <- opls(sacurine$dataMatrix, genderFc)

## ----oplsda-------------------------------------------------------------------
sacurine.oplsda <- opls(sacurine$dataMatrix, genderFc,
                        predI = 1, orthoI = NA)

## ----oplsda_subset, warning=FALSE---------------------------------------------
sacurine.oplsda <- opls(sacurine$dataMatrix, genderFc,
                        predI = 1, orthoI = NA,
                        subset = "odd")

## ----train--------------------------------------------------------------------
trainVi <- getSubsetVi(sacurine.oplsda)
confusion_train.tb <- table(genderFc[trainVi], fitted(sacurine.oplsda))
confusion_train.tb

## ----test---------------------------------------------------------------------
confusion_test.tb <- table(genderFc[-trainVi],
                           predict(sacurine.oplsda,
                                   sacurine$dataMatrix[-trainVi, ]))
confusion_test.tb

## ----get_se-------------------------------------------------------------------
data(sacurine)
sac.se <- sacurine$se

## ----se_plsda, echo = TRUE, results = "hide"----------------------------------
sac.se <- opls(sac.se, "gender")

## ----se_updated, message = FALSE----------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
message("colData:\n")
head(SummarizedExperiment::colData(sac.se))
message("\nrowData:\n")
head(SummarizedExperiment::rowData(sac.se))

## ----se_model-----------------------------------------------------------------
sac_opls.ls <- getOpls(sac.se)
names(sac_opls.ls)
plot(sac_opls.ls[["gender_PLSDA"]], typeVc = "x-score")

## ----plot_score_se, message=FALSE, warning=FALSE------------------------------
plot_score(sac.se, model.c = "gender_PLSDA", plotly.l = TRUE, info.vc = "all")

## ----get_eset-----------------------------------------------------------------
data("sacurine")
sac.set <- sacurine$eset
# viewing the ExpressionSet
# ropls::view(sac.set)

## ----eset_plsda, echo = TRUE, results = "hide"--------------------------------
# performing the PLS-DA
sac.plsda <- opls(sac.set, "gender")

## ----eset_updated-------------------------------------------------------------
sac.set <- getEset(sac.plsda)
library(Biobase)
head(Biobase::pData(sac.set))

## ----nci60_mae, message = FALSE-----------------------------------------------
data("NCI60")
nci.mae <- NCI60[["mae"]]

## ----mae_plsda----------------------------------------------------------------
nci.mae <- opls(nci.mae, "cancer",
                predI = 2, fig.pdfC = "none")

## ----mae_updated--------------------------------------------------------------
SummarizedExperiment::colData(nci.mae[["agilent"]])

## ----mae_model----------------------------------------------------------------
mae_opls.ls <- getOpls(nci.mae)
names(mae_opls.ls)
plot(mae_opls.ls[["agilent"]][["cancer_PLSDA"]], typeVc = "x-score")

## ----get_mds------------------------------------------------------------------
data("NCI60")
nci.mds <- NCI60[["mds"]]

## ----mds_plsda, echo = TRUE, results = "hide"---------------------------------
# Restricting to the "agilent" and "hgu95" datasets
nci.mds <- nci.mds[, c("agilent", "hgu95")]
# Restricting to the 'ME' and 'LE' cancer types
library(Biobase)
sample_names.vc <- Biobase::sampleNames(nci.mds[["agilent"]])
cancer_type.vc <- Biobase::pData(nci.mds[["agilent"]])[, "cancer"]
nci.mds <- nci.mds[sample_names.vc[cancer_type.vc %in% c("ME", "LE")], ]
# Building PLS-DA models for the cancer type
nci.plsda <- ropls::opls(nci.mds, "cancer", predI = 2)

## ----mds_getmset--------------------------------------------------------------
nci.mds <- ropls::getMset(nci.plsda)

## ----phenomis, eval = FALSE---------------------------------------------------
#  library(phenomis)
#  sacurine.se <- sacurine$se
#  phenomis::writing(sacurine.se, dir.c = getwd())

## ----overfit, echo = FALSE----------------------------------------------------
set.seed(123)
obsI <- 20
featVi <- c(2, 20, 200)
featMaxI <- max(featVi)
xRandMN <- matrix(runif(obsI * featMaxI), nrow = obsI)
yRandVn <- sample(c(rep(0, obsI / 2), rep(1, obsI / 2)))

layout(matrix(1:4, nrow = 2, byrow = TRUE))
for (featI in featVi) {
  randPlsi <- opls(xRandMN[, 1:featI], yRandVn,
                   predI = 2,
                   permI = ifelse(featI == featMaxI, 100, 0),
                   fig.pdfC = "none",
                   info.txtC = "none")
  plot(randPlsi, typeVc = "x-score",
       parCexN = 1.3, parTitleL = FALSE,
       parCexMetricN = 0.5)
  mtext(featI/obsI, font = 2, line = 2)
  if (featI == featMaxI)
    plot(randPlsi,
         typeVc = "permutation",
         parCexN = 1.3)
}
mtext(" obs./feat. ratio:",
      adj = 0, at = 0, font = 2,
      line = -2, outer = TRUE)

## ----vip----------------------------------------------------------------------
ageVn <- sacurine$sampleMetadata[, "age"]

pvaVn <- apply(sacurine$dataMatrix, 2,
               function(feaVn) cor.test(ageVn, feaVn)[["p.value"]])

vipVn <- getVipVn(opls(sacurine$dataMatrix, ageVn,
                       predI = 1, orthoI = NA,
                       fig.pdfC = "none"))

quantVn <- qnorm(1 - pvaVn / 2)
rmsQuantN <- sqrt(mean(quantVn^2))

opar <- par(font = 2, font.axis = 2, font.lab = 2,
            las = 1,
            mar = c(5.1, 4.6, 4.1, 2.1),
            lwd = 2, pch = 16)

plot(pvaVn, vipVn,
     col = "red",
     pch = 16,
     xlab = "p-value", ylab = "VIP", xaxs = "i", yaxs = "i")

box(lwd = 2)

curve(qnorm(1 - x / 2) / rmsQuantN, 0, 1, add = TRUE, col = "red", lwd = 3)

abline(h = 1, col = "blue")
abline(v = 0.05, col = "blue")

par(opar)

## ----se_build-----------------------------------------------------------------
# Preparing the data (matrix) and sample and variable metadata (data frames):
data(sacurine, package = "ropls")
data.mn <- sacurine$dataMatrix # matrix: samples x variables
samp.df <- sacurine$sampleMetadata # data frame: samples x sample metadata
feat.df <- sacurine$variableMetadata # data frame: features x feature metadata

# Creating the SummarizedExperiment (package SummarizedExperiment)
library(SummarizedExperiment)
sac.se <- SummarizedExperiment(assays = list(sacurine = t(data.mn)),
                               colData = samp.df,
                               rowData = feat.df)
# note that colData and rowData main format is DataFrame, but data frames are accepted when building the object
stopifnot(validObject(sac.se))

# Viewing the SummarizedExperiment
# ropls::view(sac.se)

## ----mae_build_load-----------------------------------------------------------
data("NCI60_4arrays", package = "omicade4")

## ----mae_build, message = FALSE, warning=FALSE--------------------------------
library(MultiAssayExperiment)
# Building the individual SummarizedExperiment instances
experiment.ls <- list()
sampleMap.ls <- list()
for (set.c in names(NCI60_4arrays)) {
  # Getting the data and metadata
  assay.mn <- as.matrix(NCI60_4arrays[[set.c]])
  coldata.df <- data.frame(row.names = colnames(assay.mn),
                           .id = colnames(assay.mn),
                           stringsAsFactors = FALSE) # the 'cancer' information will be stored in the main colData of the mae, not the individual SummarizedExperiments
  rowdata.df <- data.frame(row.names = rownames(assay.mn),
                           name = rownames(assay.mn),
                           stringsAsFactors = FALSE)
  # Building the SummarizedExperiment
  assay.ls <- list(se = assay.mn)
  names(assay.ls) <- set.c
  se <- SummarizedExperiment(assays = assay.ls,
                             colData = coldata.df,
                             rowData = rowdata.df)
  experiment.ls[[set.c]] <- se
  sampleMap.ls[[set.c]] <- data.frame(primary = colnames(se),
                                      colname = colnames(se)) # both datasets use identical sample names
}
sampleMap <- listToMap(sampleMap.ls)

# The sample metadata are stored in the colData data frame (both datasets have the same samples)
stopifnot(identical(colnames(NCI60_4arrays[[1]]),
                    colnames(NCI60_4arrays[[2]])))
sample_names.vc <- colnames(NCI60_4arrays[[1]])
colData.df <- data.frame(row.names = sample_names.vc,
                         cancer = substr(sample_names.vc, 1, 2))

nci.mae <- MultiAssayExperiment(experiments = experiment.ls,
                                colData = colData.df,
                                sampleMap = sampleMap)

stopifnot(validObject(nci.mae))

## ----eset_build, message = FALSE, warning = FALSE-----------------------------
# Preparing the data (matrix) and sample and variable metadata (data frames):
data(sacurine)
data.mn <- sacurine$dataMatrix # matrix: samples x variables
samp.df <- sacurine$sampleMetadata # data frame: samples x sample metadata
feat.df <- sacurine$variableMetadata # data frame: features x feature metadata
# Creating the ExpressionSet (package Biobase)
sac.set <- Biobase::ExpressionSet(assayData = t(data.mn))
Biobase::pData(sac.set) <- samp.df
Biobase::fData(sac.set) <- feat.df
stopifnot(validObject(sac.set))
# Viewing the ExpressionSet
# ropls::view(sac.set)

## ----mset_build_load----------------------------------------------------------
data("NCI60_4arrays", package = "omicade4")

## ----mset_build, message = FALSE, warning=FALSE-------------------------------
library(MultiDataSet)
# Creating the MultiDataSet instance
nci.mds <- MultiDataSet::createMultiDataSet()
# Adding the two datasets as ExpressionSet instances
for (set.c in names(NCI60_4arrays)) {
  # Getting the data
  expr.mn <- as.matrix(NCI60_4arrays[[set.c]])
  pdata.df <- data.frame(row.names = colnames(expr.mn),
                        cancer = substr(colnames(expr.mn), 1, 2),
                        stringsAsFactors = FALSE)
  fdata.df <- data.frame(row.names = rownames(expr.mn),
                        name = rownames(expr.mn),
                        stringsAsFactors = FALSE)
  # Building the ExpressionSet
  eset <- Biobase::ExpressionSet(assayData = expr.mn,
                                 phenoData = new("AnnotatedDataFrame",
                                                 data = pdata.df),
                                 featureData = new("AnnotatedDataFrame",
                                                   data = fdata.df),
                                 experimentData = new("MIAME",
                                                      title = set.c))
  # Adding to the MultiDataSet
  nci.mds <- MultiDataSet::add_eset(nci.mds, eset, dataset.type = set.c,
                                     GRanges = NA, warnings = FALSE)
}
stopifnot(validObject(nci.mds))

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

