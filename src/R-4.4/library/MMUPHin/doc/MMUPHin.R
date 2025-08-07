## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache = FALSE)

## ----Installation, eval = FALSE-----------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("MMUPHin")

## ----message=FALSE------------------------------------------------------------
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)

## ----load data----------------------------------------------------------------
data("CRC_abd", "CRC_meta")
CRC_abd[1:5, 1, drop = FALSE]

CRC_meta[1, 1:5]

table(CRC_meta$studyID)


## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics('../man/figures/adj_batch_example.png')

## ----adjust_batch-------------------------------------------------------------

fit_adjust_batch <- adjust_batch(feature_abd = CRC_abd,
                                 batch = "studyID",
                                 covariates = "study_condition",
                                 data = CRC_meta,
                                 control = list(verbose = FALSE))

CRC_abd_adj <- fit_adjust_batch$feature_abd_adj

## ----permanova----------------------------------------------------------------
library(vegan, quietly = TRUE)

D_before <- vegdist(t(CRC_abd))
D_after <- vegdist(t(CRC_abd_adj))

set.seed(1)
fit_adonis_before <- adonis2(D_before ~ studyID, data = CRC_meta)
fit_adonis_after <- adonis2(D_after ~ studyID, data = CRC_meta)
print(fit_adonis_before)
print(fit_adonis_after)

## ----lm_meta------------------------------------------------------------------
fit_lm_meta <- lm_meta(feature_abd = CRC_abd_adj,
                       batch = "studyID",
                       exposure = "study_condition",
                       covariates = c("gender", "age", "BMI"),
                       data = CRC_meta,
                       control = list(verbose = FALSE))

meta_fits <- fit_lm_meta$meta_fits

## ----significant differential abundant species--------------------------------
meta_fits %>% 
  filter(qval.fdr < 0.05) %>% 
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature)) +
  geom_bar(stat = "identity") +
  coord_flip()

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics('../man/figures/structure_types.PNG')

## ----discrete_discover--------------------------------------------------------

control_meta <- subset(CRC_meta, study_condition == "control")
control_abd_adj <- CRC_abd_adj[, rownames(control_meta)]

D_control <- vegdist(t(control_abd_adj))
fit_discrete <- discrete_discover(D = D_control,
                                  batch = "studyID",
                                  data = control_meta,
                                  control = list(k_max = 8,
                                                 verbose = FALSE))

## ----visualize discrete structure---------------------------------------------
study_id = "ZellerG_2014.metaphlan_bugs_list.stool"

internal <- data.frame(K = 2:8,
                       statistic = fit_discrete$internal_mean[, study_id],
                       se = fit_discrete$internal_se[, study_id],
                       type = "internal")

external <- data.frame(K = 2:8,
                       statistic = fit_discrete$external_mean[, study_id],
                       se = fit_discrete$external_se[, study_id],
                       type = "external")

rbind(internal, external) %>% 
  ggplot(aes(x = K, y = statistic, color = type)) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_line(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                position = position_dodge(width = 0.5), width = 0.5) +
  ggtitle("Evaluation of discrete structure in control stool microbiome (ZellerG_2014)")

## ----discrete_discover vaginal------------------------------------------------
data("vaginal_abd", "vaginal_meta")
D_vaginal <- vegdist(t(vaginal_abd))

fit_discrete_vag <- discrete_discover(D = D_vaginal,
                                      batch = "studyID",
                                      data = vaginal_meta,
                                      control = list(verbose = FALSE,
                                                     k_max = 8))

hmp_id = "HMP_2012.metaphlan_bugs_list.vagina"
data.frame(K = 2:8,
           statistic = fit_discrete_vag$internal_mean[, hmp_id],
           se = fit_discrete_vag$internal_se[, hmp_id]) %>% 
  ggplot(aes(x = K, y = statistic)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = statistic - se, 
                    ymax = statistic + se), 
                width = 0.5) +
  ggtitle("Evaluation of discrete structure in vaginal microbiome (HMP_2012)")

## ----continuos_structre-------------------------------------------------------
fit_continuous <- continuous_discover(feature_abd = control_abd_adj,
                                      batch = "studyID",
                                      data = control_meta,
                                      control = list(var_perc_cutoff = 0.5,
                                                     verbose = FALSE))

## ----visualize continuous structure-------------------------------------------
loading <- data.frame(feature = rownames(fit_continuous$consensus_loadings),
                      loading1 = fit_continuous$consensus_loadings[, 1])

loading %>%
    arrange(-abs(loading1)) %>%
    slice(1:20) %>%
    arrange(loading1) %>%
    mutate(feature = factor(feature, levels = feature)) %>%
    ggplot(aes(x = feature, y = loading1)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    ggtitle("Features with top loadings")

mds <- cmdscale(d = D_control)
colnames(mds) <- c("Axis1", "Axis2")
as.data.frame(mds) %>% 
  mutate(score1 = fit_continuous$consensus_scores[, 1]) %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = score1)) +
  geom_point() +
  coord_fixed()

## ----sessioninfo--------------------------------------------------------------
sessionInfo()

