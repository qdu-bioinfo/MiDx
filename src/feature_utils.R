# ============================
# File:feature_utils.R
# ============================
################################################################################
# Microbiome Data Preprocessing and Analysis Script
################################################################################
# Required libraries
required_packages <- c("dplyr", "ggplot2", "Maaslin2", "microeco", "magrittr", "pacman", "tidyverse","SIAMCAT","phyloseq","sva","R.utils","ANCOMBC")
lapply(required_packages, function(pkg) {
  library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
})

get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getSourceEditorContext()$path
    return(dirname(script_path))
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- sub("--file=", "", file_arg[1])
      return(dirname(script_path))
    } else {
      warning("Cannot determine script path. Returning working directory.")
      return(getwd())
    }
  }
}


################################################################################
# Main Function: analyze_microbiome_data
################################################################################
analyze_microbiome_data <- function(meta.all.filtered, feat.all, output_dir) {
 run_times <- list()
  # Function to create output folders
  create_output_folder <- function(path) {
    if (!dir.exists(path)) {
      success <- dir.create(path, recursive = TRUE)
      if (!success) {
        stop(sprintf("Unable to create directory: %s", path))
      }
    }
  }
  # Create output directories
  create_output_folder(output_dir)
  create_output_folder(file.path(output_dir, "Lefse"))
  create_output_folder(file.path(output_dir, "ANCOMBC"))
  create_output_folder(file.path(output_dir, "Maaslin2"))
  create_output_folder(file.path(output_dir, "Maaslin2", "figures"))
  # Perform Maaslin2 analysis
  start_time <- Sys.time()
  fit_data <- Maaslin2(
    input_data = feat.all,
    input_metadata = meta.all.filtered,
    transform = "none",
    output = file.path(output_dir, "Maaslin2"),
    fixed_effects = c("Group"),
    correction = "BH"
  )
  end_time <- Sys.time()
  run_times$Maaslin2 <- as.numeric(difftime(end_time, start_time, units = "secs"))
  # Lefse analysis
  tax <- data.frame(
    kingdom = NA,
    Phylum = NA,
    Class = NA,
    Order = NA,
    Family = NA,
    Genus = NA,
    Species = row.names(feat.all)
  )
  row.names(tax) <- row.names(feat.all)
  dataset <- microtable$new(
    sample_table = meta.all.filtered,
    otu_table = as.data.frame(feat.all),
    tax_table = tax
  )
  start_time <- Sys.time()
  lefse <- trans_diff$new(
    dataset = dataset,
    method = "lefse",
    group = "Group",
    alpha = 1,
    taxa_level = "species",
    p_adjust_method = "fdr"
  )
  end_time <- Sys.time()
  run_times$Lefse <- as.numeric(difftime(end_time, start_time, units = "secs"))
  # Save Lefse results
  write.table(
    lefse$res_diff,
    file.path(output_dir, "Lefse/lefse.csv"),
    sep = ','
  )
  set.seed(42)

    otu <- otu_table(as.matrix(feat.all), taxa_are_rows = TRUE)

    meta <- sample_data(meta.all.filtered)

    ps <- phyloseq(otu, meta)

  start_time <- Sys.time()
  out <- ancombc2(
    data           = ps,
    #taxa_are_rows  = TRUE,
    #meta_data      = meta.all.filtered,
    fix_formula    = "Group",
    rand_formula   = NULL,
    p_adj_method   = "BH",
    prv_cut        = 0.10,
    lib_cut        = 0,
    s0_perc        = 0.05,
    group          = "Group",
    struc_zero     = TRUE,
    neg_lb         = TRUE,
    alpha          = 0.05,
    n_cl           = 1,
    verbose        = TRUE,
    global         = FALSE,
    pairwise       = FALSE,
    dunnet         = FALSE,
    trend          = FALSE
  )
  end_time <- Sys.time()
  run_times$ANCOMBC <- as.numeric(difftime(end_time, start_time, units = "secs"))
  result_df <- data.frame(
    taxon = out$res$taxon,
    W     = out$res$W_GroupCTR,
    pval  = out$res$p_GroupCTR,
    qval  = out$res$q_GroupCTR,
    sig   = out$res$diff_GroupCTR
  )
  write.csv(
    result_df,
    file.path(output_dir, "ANCOMBC/ancombc.csv"),
    row.names = FALSE,quote     = FALSE
  )
  # Function to run external analysis tools
  run_analysis_tools <- function(feat.all, meta.all.filtered, output_dir) {

    script_dir <- get_script_dir()
    print(script_dir)
    tool_script_dir <- file.path(script_dir, "Tool_script")
    analysis_tools <- list(
    "Wilcoxon" = file.path(tool_script_dir,"Run_Wilcox.R") ,
    "metagenomeSeq" = file.path(tool_script_dir,"Run_metagenomeSeq.R"),
    "t_test" = file.path(tool_script_dir,"Run_t_test.R")
    #"ANCOM" = paste(script_dir, "/Tool_script/Run_ANCOM.R", sep = "")
    )
    for (method in names(analysis_tools)) {
      tool_script <- analysis_tools[[method]]
      output_path <- normalizePath(file.path(output_dir, method, paste(method, ".tsv", sep = "")), mustWork = FALSE)
      dir_path <- dirname(output_path)
      # Create directories if needed
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
        cat("Directory created: ", dir_path, "\n")
      }
        # ---- in run_analysis_tools() ----

        tmp_dir <- file.path(output_dir, "temporary")
        if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)


        feat_file <- file.path(tmp_dir, "feat.csv")
        meta_file <- file.path(tmp_dir, "meta.csv")
        write.csv(feat.all, feat_file, row.names = TRUE)
        write.csv(meta.all.filtered, meta_file, row.names = TRUE)

      tryCatch({
        if (method == "ANCOM") {
          ANCOM_DIR <- paste("./Tool_script/Ancom2_Script", sep = "")
          ANCOM_Script_Path <- file.path(ANCOM_DIR, "ancom_v2.1.R")
          cmd <- sprintf("Rscript \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"", tool_script, feat_file, meta_file, output_path, ANCOM_Script_Path)
        } else {
          cmd <- sprintf("Rscript \"%s\" \"%s\" \"%s\" \"%s\"", tool_script, feat_file, meta_file, output_path)
        }
        # Execute command
        start_time <- Sys.time()
        system(cmd)
        end_time <- Sys.time()
        run_times[[method]] <<- as.numeric(difftime(end_time, start_time, units = "secs"))
        cat(sprintf("Method %s executed successfully.\n", method))
      }, error = function(e) {
        cat(sprintf("Error executing method %s: %s\n", method, e$message))
        stop(sprintf("Execution halted: Method %s failed.", method))
      })
    }
  }
  # Run external analysis tools
  run_analysis_tools(feat.all, meta.all.filtered, output_dir)
  write.csv(
  data.frame(method = names(run_times), seconds = unlist(run_times)),
  file = file.path(output_dir, "runtime_summary.csv"),
  row.names = FALSE
  )
  sprintf(file.path(output_dir, "runtime_summary.csv"))
}

read_table_and_check_line_count <- function(filepath, ...) {

  raw_lines <- readLines(filepath, warn = FALSE)
  lines <- raw_lines[ nzchar(raw_lines) ]
  total_lines <- length(lines)

  header_present <- {
    first_row <- lines[1]

    cols <- strsplit(first_row, "[,\t]")[[1]]
    !all(grepl("^V\\d+$", cols))
  }
  expected_rows <- if (header_present) total_lines - 1L else total_lines

  df <- read.table(filepath, ...)

  if (nrow(df) != expected_rows) {
    stop(sprintf(
      "Expected %d data rows (from %d non-blank lines), but found %d rows.",
      expected_rows, total_lines, nrow(df)
    ))
  }
  df

}

read_genera_hackathon_results <- function(feature_dir) {

  da_tool_filepath <- list()
  #da_tool_filepath[["ancom"]] <- paste(feature_dir,"/ANCOM/","Ancom.tsv", sep = "")
  da_tool_filepath[["maaslin2"]] <- paste(feature_dir,"/Maaslin2/all_results.tsv", sep = "")
  da_tool_filepath[["metagenomeSeq"]] <- paste(feature_dir,"/metagenomeSeq/metagenomeSeq.tsv", sep = "")
  da_tool_filepath[["ttest"]] <- paste(feature_dir,"/t_test/t_test.tsv", sep = "")
  da_tool_filepath[["wilcoxon"]] <- paste(feature_dir,"/Wilcoxon/Wilcoxon.tsv", sep = "")
  da_tool_filepath[["lefse"]] <- paste(feature_dir,"/Lefse/lefse.csv", sep = "")
  da_tool_filepath[["ancombc"]] <- paste(feature_dir,"/ANCOMBC/ancombc.csv", sep = "")

  adjP_colname <- list()
  #adjP_colname[["ancom"]] <- "detected_0.8"
  adjP_colname[["maaslin2"]] <- "qval"
  adjP_colname[["metagenomeSeq"]] <- "adjPvalues"
  adjP_colname[["ttest"]] <- "x"
  adjP_colname[["wilcoxon"]] <- "x"
  adjP_colname[["lefse"]] <- "P.adj"
  adjP_colname[["ancombc"]] <- "qval"

  da_tool_results <- list()
  missing_tools <- c()
  for(da_tool in names(da_tool_filepath)) {
    if(! (file.exists(da_tool_filepath[[da_tool]]))) {
      missing_tools <- c(missing_tools, da_tool)
      message(paste("File ", da_tool_filepath[[da_tool]], " not found. Skipping.", sep=""))
      next
    }
    if(da_tool %in% c("ancombc")) {
      da_tool_results[["ancombc"]] <- read_table_and_check_line_count(da_tool_filepath[[da_tool]], sep=",", row.names=1, header=TRUE,fill=TRUE,quote="")
    } else if(da_tool%in% c("lefse")) {
      da_tool_results[[da_tool]] <- read_table_and_check_line_count(da_tool_filepath[[da_tool]], sep=",", row.names=1,header=TRUE, stringsAsFactors=FALSE)
      da_tool_results[[da_tool]] <- da_tool_results[[da_tool]][da_tool_results[[da_tool]][["LDA"]] >= 2, , drop = FALSE]
      row_names <- rownames(da_tool_results[[da_tool]])
      row_names <- sub(".*\\|", "", row_names)
      rownames(da_tool_results[[da_tool]]) <- row_names
    }else if(da_tool%in% c("wilcoxon","ttest")) {
      da_tool_results[[da_tool]] <- read_table_and_check_line_count(da_tool_filepath[[da_tool]], sep="\t", row.names=1, header=TRUE,fill=TRUE,quote="")
    } else {
      da_tool_results[[da_tool]] <- read_table_and_check_line_count(da_tool_filepath[[da_tool]], sep="\t", row.names=1, header=TRUE,fill=TRUE,quote="")
    }
  }

  for(da_tool in names(da_tool_filepath)) {
    rownames(da_tool_results[[da_tool]]) <- gsub("\\.", "_", rownames(da_tool_results[[da_tool]]))
    rownames(da_tool_results[[da_tool]]) <- gsub("/", "_", rownames(da_tool_results[[da_tool]]))
  }

  all_rows <- c()
  for(da_tool in names(adjP_colname)) {
    all_rows <- c(all_rows, rownames(da_tool_results[[da_tool]]))
  }
  all_rows <- all_rows[-which(duplicated(all_rows))]
  adjP_table <- data.frame(matrix(NA, ncol=length(names(da_tool_results)), nrow=length(all_rows)))
  colnames(adjP_table) <- names(da_tool_results)
  rownames(adjP_table) <- all_rows
  for(da_tool in colnames(adjP_table)) {
    if(da_tool %in% missing_tools) {
      next
    }
    if(da_tool == "lefse") {
      adjP_table[rownames(da_tool_results[[da_tool]]), da_tool] <-da_tool_results[[da_tool]][, adjP_colname[[da_tool]]]
    } else if(da_tool =="ancom") {
      sig_ancom_hits <- which(da_tool_results[[da_tool]][, adjP_colname[[da_tool]]])
      ancom_results <- rep(1, length(da_tool_results[[da_tool]][, adjP_colname[[da_tool]]]))
      ancom_results[sig_ancom_hits] <- 0
      adjP_table[rownames(da_tool_results[[da_tool]]), da_tool] <- ancom_results
    }else if (da_tool%in% c("wilcoxon","ttest")){
      adjP_table[rownames(da_tool_results[[da_tool]]), da_tool] <- p.adjust(da_tool_results[[da_tool]][, adjP_colname[[da_tool]]], "fdr")
    }else {
      adjP_table[rownames(da_tool_results[[da_tool]]), da_tool] <- da_tool_results[[da_tool]][, adjP_colname[[da_tool]]]
    }
  }
  return(list(raw_tables=da_tool_results,adjP_table=adjP_table))
}
    run_feature_analysis <-  function(
                                 norm_mat,
                                 meta_filtered,
                                 feature_script_dir

    ) {
    output_dir <- file.path(feature_script_dir,"All_features_tools")
    analyze_microbiome_data(meta_filtered,norm_mat,output_dir)

    merge_feature <- read_genera_hackathon_results(output_dir)
    raw_tables<-merge_feature$raw_tables
    adjP_table<-merge_feature$adjP_table
    return (adjP_table)
    }
