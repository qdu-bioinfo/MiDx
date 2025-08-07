## ----load_files, message=FALSE------------------------------------------------
library("SIAMCAT")

data("feat_crc_zeller", package="SIAMCAT")
data("meta_crc_zeller", package="SIAMCAT")

## ----show_features------------------------------------------------------------
feat.crc.zeller[1:3, 1:3]
dim(feat.crc.zeller)

## ----show_meta----------------------------------------------------------------
head(meta.crc.zeller)

## ----create_label-------------------------------------------------------------
label.crc.zeller <- create.label(meta=meta.crc.zeller,
    label='Group', case='CRC')

## ----start--------------------------------------------------------------------
sc.obj <- siamcat(feat=feat.crc.zeller,
    label=label.crc.zeller,
    meta=meta.crc.zeller)

## ----show_siamcat-------------------------------------------------------------
show(sc.obj)

## ----filter_feat--------------------------------------------------------------
sc.obj <- filter.features(sc.obj,
    filter.method = 'abundance',
    cutoff = 0.001)

## ----check_associations, eval=FALSE-------------------------------------------
#  sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
#  association.plot(sc.obj, sort.by = 'fc',
#                  panels = c('fc', 'prevalence', 'auroc'))

## ----check_confounders, eval=FALSE--------------------------------------------
#  check.confounders(sc.obj, fn.plot = 'confounder_plots.pdf',
#                      meta.in = NULL, feature.type = 'filtered')

## ----normalize_feat-----------------------------------------------------------
sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
    norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))

## ----data_split---------------------------------------------------------------
sc.obj <-  create.data.split(sc.obj, num.folds = 5, num.resample = 2)

## ----train_model, message=FALSE, results='hide'-------------------------------
sc.obj <- train.model(sc.obj, method = "lasso")

## ----show_models--------------------------------------------------------------
# get information about the model type
model_type(sc.obj)

# access the models
models <- models(sc.obj)
models[[1]]$model

## ----make_predictions, message=FALSE, results='hide'--------------------------
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)

## ----pred_matrix_head---------------------------------------------------------
head(pred_matrix)

## ----eval_predictions---------------------------------------------------------
sc.obj <-  evaluate.predictions(sc.obj)

## ----eval_plot, eval=FALSE----------------------------------------------------
#  model.evaluation.plot(sc.obj)

## ----eval=FALSE---------------------------------------------------------------
#  model.interpretation.plot(sc.obj, fn.plot = 'interpretation.pdf',
#      consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

## -----------------------------------------------------------------------------
sessionInfo()

