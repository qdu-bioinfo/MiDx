## ----setup, message = FALSE, warning = FALSE, comment = NA--------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
                      fig.width = 6.25, fig.height = 5)
library(ANCOMBC)
library(tidyverse)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
  "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

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
out = ancombc(data = pseq, tax_level = "Family", 
              formula = "age + region + bmi", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)

res = out$res
res_global = out$res_global

## ----eval=FALSE---------------------------------------------------------------
# out = ancombc(data = pseq, tax_level = "Family",
#               formula = "age + region * bmi",
#               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#               group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
#               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
#               n_cl = 1, verbose = TRUE)
# 
# res = out$res
# res_global = out$res_global

## -----------------------------------------------------------------------------
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "Age", "NE - CE", "SE - CE", 
             "US - CE", "Overweight - Obese", "Lean - Obese")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

## -----------------------------------------------------------------------------
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

## -----------------------------------------------------------------------------
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

## -----------------------------------------------------------------------------
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

## -----------------------------------------------------------------------------
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

## -----------------------------------------------------------------------------
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

## -----------------------------------------------------------------------------
samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log observed abundances")

## -----------------------------------------------------------------------------
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
    mutate(taxon_id = res$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
    dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_age = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, age, ageSE) %>%
  dplyr::filter(age != 0) %>% 
  dplyr::arrange(desc(age)) %>%
  dplyr::mutate(direct = ifelse(age > 0, "Positive LFC", "Negative LFC"))
df_fig_age$taxon_id = factor(df_fig_age$taxon_id, levels = df_fig_age$taxon_id)
df_fig_age$direct = factor(df_fig_age$direct, 
                        levels = c("Positive LFC", "Negative LFC"))
  
p_age = ggplot(data = df_fig_age, 
           aes(x = taxon_id, y = age, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = age - ageSE, ymax = age + ageSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_age

## -----------------------------------------------------------------------------
df_fig_bmi = df_lfc %>% 
  filter(bmioverweight != 0 | bmilean != 0) %>%
  transmute(taxon_id, 
            `Overweight vs. Obese` = round(bmioverweight, 2),
            `Lean vs. Obese` = round(bmilean, 2)) %>%
  pivot_longer(cols = `Overweight vs. Obese`:`Lean vs. Obese`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_bmi$value))
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2
p_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_bmi

## -----------------------------------------------------------------------------
tab_w = res_global[, c("taxon", "W")]
tab_w %>% datatable(caption = "Test Statistics 
                    from the Global Test Result") %>%
      formatRound(c("W"), digits = 2)

## -----------------------------------------------------------------------------
tab_p = res_global[, c("taxon", "p_val")]
tab_p %>% datatable(caption = "P-values 
                    from the Global Test Result") %>%
      formatRound(c("p_val"), digits = 2)

## -----------------------------------------------------------------------------
tab_q = res_global[, c("taxon", "q_val")]
tab_q %>% datatable(caption = "Adjusted p-values 
                    from the Global Test Result") %>%
      formatRound(c("q_val"), digits = 2)

## -----------------------------------------------------------------------------
tab_diff = res_global[, c("taxon", "diff_abn")]
tab_diff %>% datatable(caption = "Differentially Abundant Taxa 
                       from the Global Test Result")

## -----------------------------------------------------------------------------
sig_taxa = res_global %>%
  dplyr::filter(diff_abn == TRUE) %>%
  .$taxon

df_bmi = tab_lfc %>%
    dplyr::select(Taxon, `Overweight - Obese`, `Lean - Obese`) %>%
    filter(Taxon %in% sig_taxa)

df_heat = df_bmi %>%
    pivot_longer(cols = -one_of("Taxon"),
                 names_to = "region", values_to = "value") %>%
    mutate(value = round(value, 2))
df_heat$Taxon = factor(df_heat$Taxon, levels = sort(sig_taxa))

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
p_heat = df_heat %>%
  ggplot(aes(x = region, y = Taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, Taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, 
       title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_heat

## ----eval=FALSE---------------------------------------------------------------
# tse = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)
# 
# out = ancombc(data = tse, assay_name = "counts", tax_level = "Family",
#               formula = "age + region + bmi",
#               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#               group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
#               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
#               n_cl = 1, verbose = TRUE)
# 
# res = out$res
# res_global = out$res_global

## ----eval=FALSE---------------------------------------------------------------
# abundance_data = microbiome::abundances(pseq)
# aggregate_data = microbiome::abundances(microbiome::aggregate_taxa(pseq, "Family"))
# meta_data = microbiome::meta(pseq)
# 
# out = ancombc(data = abundance_data, aggregate_data = aggregate_data,
#               meta_data = meta_data, formula = "age + region + bmi",
#               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#               group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
#               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
#               n_cl = 1, verbose = TRUE)
# 
# res = out$res
# res_global = out$res_global

## ----sessionInfo, message = FALSE, warning = FALSE, comment = NA--------------
sessionInfo()

