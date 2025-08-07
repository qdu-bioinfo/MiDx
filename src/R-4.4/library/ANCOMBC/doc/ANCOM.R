## ----setup, message = FALSE, warning = FALSE, comment = NA--------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA,
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

# Subset to baseline
pseq = phyloseq::subset_samples(atlas1006, time == 0)

# Re-code the bmi group
meta_data = microbiome::meta(pseq)
meta_data$bmi = recode(meta_data$bmi_group,
                       obese = "obese",
                       severeobese = "obese",
                       morbidobese = "obese")

# Note that by default, levels of a categorical variable in R are sorted 
# alphabetically. In this case, the reference level for `bmi` will be 
# `lean`. To manually change the reference level, for instance, setting `obese`
# as the reference level, use:
meta_data$bmi = factor(meta_data$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(meta_data$bmi)

# Create the region variable
meta_data$region = recode(as.character(meta_data$nationality),
                          Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                          CentralEurope = "CE", EasternEurope = "EE",
                          .missing = "unknown")

phyloseq::sample_data(pseq) = meta_data

# Subset to lean, overweight, and obese subjects
pseq = phyloseq::subset_samples(pseq, bmi %in% c("lean", "overweight", "obese"))
# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
pseq = phyloseq::subset_samples(pseq, ! region %in% c("EE", "unknown"))

print(pseq)

## -----------------------------------------------------------------------------
set.seed(123)
out = ancom(data = pseq, tax_level = "Family", meta_data = NULL,
            p_adj_method = "holm", prv_cut = 0.10,
            lib_cut = 1000, main_var = "bmi", adj_formula = "age + region", 
            rand_formula = NULL, lme_control = NULL, struc_zero = TRUE,
            neg_lb = TRUE, alpha = 0.05, n_cl = 2, verbose = TRUE)

res = out$res

# Similarly, if the main variable of interest is continuous, such as age, the
# ancom model can be specified as
# out = ancom(data = pseq, tax_level = "Family", meta_data = NULL,
#             p_adj_method = "holm", prv_cut = 0.10,
#             lib_cut = 1000, main_var = "age", adj_formula = "bmi + region",
#             rand_formula = NULL, lme_control = NULL, struc_zero = FALSE,
#             neg_lb = FALSE, alpha = 0.05, n_cl = 2, verbose = TRUE)

## -----------------------------------------------------------------------------
q_val = out$q_data
beta_val = out$beta_data
# Only consider the effect sizes with the corresponding q-value less than alpha
beta_val = beta_val * (q_val < 0.05) 
# Choose the maximum of beta's as the effect size
beta_pos = apply(abs(beta_val), 2, which.max) 
beta_max = vapply(seq_along(beta_pos), function(i) 
    beta_val[beta_pos[i], i], FUN.VALUE = double(1))
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(out$zero_ind), 
                nrow(tse), 
                sum(apply(out$zero_ind[, -1], 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = 0.7 * (n_taxa - 1)

df_fig_w = res %>%
  dplyr::mutate(beta = beta_max,
                direct = case_when(
                  detected_0.7 == TRUE & beta > 0 ~ "Positive",
                  detected_0.7 == TRUE & beta <= 0 ~ "Negative",
                  TRUE ~ "Not Significant"
                  )) %>%
  dplyr::arrange(W)
df_fig_w$taxon = factor(df_fig_w$taxon, levels = df_fig_w$taxon)
df_fig_w$W = replace(df_fig_w$W, is.infinite(df_fig_w$W), n_taxa - 1)
df_fig_w$direct = factor(df_fig_w$direct, 
                         levels = c("Negative", "Positive", "Not Significant"))

p_w = df_fig_w %>%
  ggplot(aes(x = taxon, y = W, color = direct)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x = "Taxon", y = "W") +
  scale_color_discrete(name = NULL) + 
  geom_hline(yintercept = cut_off, linetype = "dotted", 
             color = "blue", size = 1.5) +
  geom_text(aes(x = 2, y = cut_off + 0.5, label = "W[0.7]"), 
            size = 5, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())
p_w

## ----eval=FALSE---------------------------------------------------------------
# tse = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
# 
# set.seed(123)
# out = ancom(data = tse, assay_name = "counts",
#             tax_level = "Family", meta_data = NULL,
#             p_adj_method = "holm", prv_cut = 0.10,
#             lib_cut = 1000, main_var = "bmi", adj_formula = "age + region",
#             rand_formula = NULL, lme_control = NULL, struc_zero = TRUE,
#             neg_lb = TRUE, alpha = 0.05, n_cl = 2, verbose = TRUE)
# 
# res = out$res

## ----eval=FALSE---------------------------------------------------------------
# abundance_data = microbiome::abundances(pseq)
# aggregate_data = microbiome::abundances(microbiome::aggregate_taxa(pseq, "Family"))
# meta_data = microbiome::meta(pseq)
# 
# set.seed(123)
# out = ancom(data = abundance_data, aggregate_data = aggregate_data,
#             meta_data = meta_data, p_adj_method = "holm", prv_cut = 0.10,
#             lib_cut = 1000, main_var = "bmi", adj_formula = "age + region",
#             rand_formula = NULL, lme_control = NULL, struc_zero = TRUE,
#             neg_lb = TRUE, alpha = 0.05, n_cl = 2, verbose = TRUE)
# 
# res = out$res

## -----------------------------------------------------------------------------
data(dietswap, package = "microbiome")
print(dietswap)

## -----------------------------------------------------------------------------
set.seed(123)
out = ancom(data = dietswap, tax_level = "Family", 
            p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
            main_var = "group",
            adj_formula = "nationality + timepoint", 
            rand_formula = "(timepoint | subject)", 
            lme_control = lme4::lmerControl(), 
            struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 2)

res = out$res

## -----------------------------------------------------------------------------
q_val = out$q_data
beta_val = out$beta_data
# Only consider the effect sizes with the corresponding q-value less than alpha
beta_val = beta_val * (q_val < 0.05) 
# Choose the maximum of beta's as the effect size
beta_pos = apply(abs(beta_val), 2, which.max) 
beta_max = vapply(seq_along(beta_pos), function(i) beta_val[beta_pos[i], i],
                  FUN.VALUE = double(1))
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(out$zero_ind), 
                nrow(tse), 
                sum(apply(out$zero_ind[, -1], 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = 0.7 * (n_taxa - 1)

df_fig_w = res %>%
  dplyr::mutate(beta = beta_max,
                direct = case_when(
                  detected_0.7 == TRUE & beta > 0 ~ "Positive",
                  detected_0.7 == TRUE & beta <= 0 ~ "Negative",
                  TRUE ~ "Not Significant"
                  )) %>%
  dplyr::arrange(W)
df_fig_w$taxon = factor(df_fig_w$taxon, levels = df_fig_w$taxon)
df_fig_w$W = replace(df_fig_w$W, is.infinite(df_fig_w$W), n_taxa - 1)
df_fig_w$direct = factor(df_fig_w$direct, 
                     levels = c("Negative", "Positive", "Not Significant"))

p_w = df_fig_w %>%
  ggplot(aes(x = taxon, y = W, color = direct)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x = "Taxon", y = "W") +
  scale_color_discrete(name = NULL) + 
  geom_hline(yintercept = cut_off, linetype = "dotted", 
             color = "blue", size = 1.5) +
  geom_text(aes(x = 2, y = cut_off + 0.5, label = "W[0.7]"), 
            size = 5, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())
p_w

## ----eval=FALSE---------------------------------------------------------------
# tse = mia::makeTreeSummarizedExperimentFromPhyloseq(dietswap)
# 
# set.seed(123)
# out = ancom(data = tse, assay_name = "counts", tax_level = "Family",
#             p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#             main_var = "group",
#             adj_formula = "nationality + timepoint",
#             rand_formula = "(timepoint | subject)",
#             lme_control = lme4::lmerControl(),
#             struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 2)
# 
# res = out$res

## ----eval=FALSE---------------------------------------------------------------
# abundance_data = microbiome::abundances(dietswap)
# aggregate_data = microbiome::abundances(microbiome::aggregate_taxa(dietswap, "Family"))
# meta_data = microbiome::meta(dietswap)
# 
# set.seed(123)
# out = ancom(data = abundance_data, aggregate_data = aggregate_data,
#             meta_data = meta_data, p_adj_method = "holm",
#             prv_cut = 0.10, lib_cut = 1000, main_var = "group",
#             adj_formula = "nationality + timepoint",
#             rand_formula = "(timepoint | subject)",
#             lme_control = lme4::lmerControl(),
#             struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 2)
# 
# res = out$res

## ----sessionInfo, message = FALSE, warning = FALSE, comment = NA--------------
sessionInfo()

