# ============================
# File: siamcat_utils.R
# ============================

library(dplyr)
library(SIAMCAT)
library(phyloseq)
library(rstudioapi)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
run_siamcat_filter_analysis <- function(raw_filt,
                                 meta_filt,
                                 groups,
                                 filter_method,
                                 thresholds_list) {
  if (!all(groups %in% unique(meta_filt$Group))) {
    stop("groups 中的某些标签在 meta_filt$Group 中未找到")
  }

  set.seed(42)
  lbl <- create.label(
    meta    = meta_filt,
    label   = "Group",
    case    = groups[2],
    control = groups[1]
  )
  sc <- siamcat(feat = t(raw_filt), label = lbl, meta = meta_filt)
    sc_f <- filter.features(
    siamcat       = sc,
    filter.method = filter_method,
    cutoff        = thresholds_list[[filter_method]],
    feature.type  = "original",
    verbose     = 2
  )
  mat <- get.filt_feat.matrix(sc_f)
  return(as.data.frame(mat))
}

run_siamcat_norm_analysis <- function(raw_filt,
                                 meta_filt,
                                 groups,
                                 norm_method) {
  if (!all(groups %in% unique(meta_filt$Group))) {
    stop("groups 中的某些标签在 meta_filt$Group 中未找到")
  }

  params <- switch(norm_method,
    log.unit = list(log.n0 = 1e-6, n.p = 2, norm.margin = 1),
    log.std  = list(log.n0 = 1e-6, sd.min.q = 0.1),
    log.clr  = list(log.n0 = 1e-6, sd.min.q = 0.1),
    std      = list(sd.min.q = 0.1),
    rank.std = list(sd.min.q = 0.1),
    list()
  )
  set.seed(42)
  lbl <- create.label(
    meta    = meta_filt,
    label   = "Group",
    case    = groups[2],
    control = groups[1]
  )
  sc <- siamcat(feat = t(raw_filt), label = lbl, meta = meta_filt)

  sc_f <- filter.features(
    siamcat       = sc,
    filter.method ="pass",
    cutoff        = 0,
    feature.type  = "original",
    verbose     = 2
  )
  sc_n <- normalize.features(
    sc_f,
    norm.method = norm_method,
    norm.param  = params,
    verbose     = 2
  )
  mat <-get.norm_feat.matrix(sc_n)
  return(as.data.frame(mat))
}

