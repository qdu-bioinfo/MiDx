## ----start, message=FALSE, warning=FALSE--------------------------------------
library("SIAMCAT")

data.loc <- 'https://zenodo.org/api/files/d81e429c-870f-44e0-a44a-2a4aa541b6c1/'
# this is data from Zeller et al., Mol. Syst. Biol. 2014
fn.meta.fr  <- paste0(data.loc, 'meta_Zeller.tsv')
fn.feat.fr  <- paste0(data.loc, 'specI_Zeller.tsv')

# this is the external dataset from Yu et al., Gut 2017
fn.meta.cn  <- paste0(data.loc, 'meta_Yu.tsv')
fn.feat.cn  <- paste0(data.loc, 'specI_Yu.tsv')

## ----siamcat_fr---------------------------------------------------------------
# features
# be vary of the defaults in R!!!
feat.fr  <- read.table(fn.feat.fr, sep='\t', quote="",
    check.names = FALSE, stringsAsFactors = FALSE)
# the features are counts, but we want to work with relative abundances
feat.fr.rel <- prop.table(as.matrix(feat.fr), 2)

# metadata
meta.fr  <- read.table(fn.meta.fr, sep='\t', quote="",
    check.names=FALSE, stringsAsFactors=FALSE)

# create SIAMCAT object
siamcat.fr <- siamcat(feat=feat.fr.rel, meta=meta.fr,
    label='Group', case='CRC')

## ----siamcat_cn---------------------------------------------------------------
# features
feat.cn  <- read.table(fn.feat.cn, sep='\t', quote="",
    check.names = FALSE)
feat.cn.rel <- prop.table(as.matrix(feat.cn), 2)

# metadata
meta.cn  <- read.table(fn.meta.cn, sep='\t', quote="",
    check.names=FALSE, stringsAsFactors = FALSE)

# SIAMCAT object
siamcat.cn <- siamcat(feat=feat.cn.rel, meta=meta.cn,
        label='Group', case='CRC')

## ----preprocessing_fr---------------------------------------------------------
siamcat.fr <- filter.features(
    siamcat.fr,
    filter.method = 'abundance',
    cutoff = 0.001,
    rm.unmapped = TRUE,
    verbose=2
)

siamcat.fr <- normalize.features(
    siamcat.fr,
    norm.method = "log.std",
    norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1),
    verbose = 2
)

## ----build_model_fr, results='hide'-------------------------------------------
siamcat.fr <-  create.data.split(
    siamcat.fr,
    num.folds = 5,
    num.resample = 2
)

siamcat.fr <- train.model(
    siamcat.fr,
    method = "lasso"
)

## ----predict_evaluate_fr, results='hide'--------------------------------------
siamcat.fr <- make.predictions(siamcat.fr)

siamcat.fr <-  evaluate.predictions(siamcat.fr)

## ----normalize_cn-------------------------------------------------------------

siamcat.cn <- normalize.features(siamcat.cn,
    norm.param=norm_params(siamcat.fr),
    feature.type='original',
    verbose = 2)


## ----predict_cn, results='hide'-----------------------------------------------
siamcat.cn <- make.predictions(
    siamcat = siamcat.fr,
    siamcat.holdout = siamcat.cn,
    normalize.holdout = FALSE)

## ----alternative_pipeline_cn, eval=FALSE--------------------------------------
#  ## Alternative Code, not run here
#  siamcat.cn <- siamcat(feat=feat.cn.rel, meta=meta.cn,
#      label='Group', case='CRC')
#  siamcat.cn <- make.predictions(siamcat = siamcat.fr,
#      siamcat.holdout = siamcat.cn,
#      normalize.holdout = TRUE)

## ----eval_cn, message=FALSE---------------------------------------------------
siamcat.cn <- evaluate.predictions(siamcat.cn)

## ----eval_plot, eval=FALSE----------------------------------------------------
#  model.evaluation.plot('FR-CRC'=siamcat.fr,
#      'CN-CRC'=siamcat.cn,
#      colours=c('dimgrey', 'orange'))

## ----session_info-------------------------------------------------------------
sessionInfo()

