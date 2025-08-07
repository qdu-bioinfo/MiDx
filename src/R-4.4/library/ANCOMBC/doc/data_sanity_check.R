## ----setup, message = FALSE, warning = FALSE, comment = NA--------------------
knitr::opts_chunk$set(warning = FALSE, comment = NA,
                      fig.width = 6.25, fig.height = 5)
library(ANCOMBC)
library(tidyverse)

## ----getPackage, eval=FALSE---------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ANCOMBC")

## ----load, eval=FALSE---------------------------------------------------------
# library(ANCOMBC)

## -----------------------------------------------------------------------------
data(atlas1006, package = "microbiome")

atlas1006

## -----------------------------------------------------------------------------
phyloseq::rank_names(atlas1006)

## -----------------------------------------------------------------------------
colnames(microbiome::meta(atlas1006))

## -----------------------------------------------------------------------------
# With `group` variable
check_results = data_sanity_check(data = atlas1006,
                                  tax_level = "Family",
                                  fix_formula = "age + sex + bmi_group",
                                  group = "bmi_group",
                                  struc_zero = TRUE,
                                  global = TRUE,
                                  verbose = TRUE)

## -----------------------------------------------------------------------------
# Without `group` variable
check_results = data_sanity_check(data = atlas1006,
                                  tax_level = "Family",
                                  fix_formula = "age + sex + bmi_group",
                                  group = NULL,
                                  struc_zero = FALSE,
                                  global = FALSE,
                                  verbose = TRUE)

## -----------------------------------------------------------------------------
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

## -----------------------------------------------------------------------------
mia::taxonomyRanks(tse)

## -----------------------------------------------------------------------------
colnames(SummarizedExperiment::colData(tse))

## -----------------------------------------------------------------------------
check_results = data_sanity_check(data = tse,
                                  assay_name = "counts",
                                  tax_level = "Family",
                                  fix_formula = "age + sex + bmi_group",
                                  group = "bmi_group",
                                  struc_zero = TRUE,
                                  global = TRUE,
                                  verbose = TRUE)

## -----------------------------------------------------------------------------
abundance_data = microbiome::abundances(atlas1006)
meta_data = microbiome::meta(atlas1006)

## -----------------------------------------------------------------------------
all(rownames(meta_data) %in% colnames(abundance_data))

## -----------------------------------------------------------------------------
colnames(meta_data)

## -----------------------------------------------------------------------------
check_results = data_sanity_check(data = abundance_data,
                                  assay_name = "counts",
                                  tax_level = "Family",
                                  meta_data = meta_data,
                                  fix_formula = "age + sex + bmi_group",
                                  group = "bmi_group",
                                  struc_zero = TRUE,
                                  global = TRUE,
                                  verbose = TRUE)

## ----sessionInfo, message = FALSE, warning = FALSE, comment = NA--------------
sessionInfo()

