# ============================
# File:correct_utils.R
# ============================

library(sva)
library(SIAMCAT)
library(phyloseq)
library(dplyr)
library(limma)
library(MMUPHin)
run_correct_analysis <- function(class_level,
                                 groups,
                                 Correct_method,
                                 norm_mat,
                                 meta_filtered
) {
    if (Correct_method == "combat") {
        expr_combat <- ComBat(dat = as.matrix(norm_mat), batch = meta_filtered$Study)
        correct_table <- otu_table(expr_combat, taxa_are_rows = TRUE)

    } else if (Correct_method == "limma") {
      expr_rbe <- removeBatchEffect(as.matrix(norm_mat), batch = meta_filtered$Study)
      correct_table <- otu_table(expr_rbe, taxa_are_rows = TRUE)

    } else if (Correct_method == "MMUPHin") {
         meta_filtered$StudyID <- factor(meta_filtered$Study)
         fit_adjust_batch <- adjust_batch(
            feature_abd = as.matrix(norm_mat),
            batch = "StudyID",
            data = meta_filtered
        )
        correct_table <- otu_table(fit_adjust_batch$feature_abd_adj, taxa_are_rows = TRUE)
    } else {
        correct_table <- norm_mat
    }
    return(as.data.frame(correct_table))
}
