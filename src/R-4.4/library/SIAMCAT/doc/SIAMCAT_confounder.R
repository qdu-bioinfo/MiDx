## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library("tidyverse")
library("SIAMCAT")
library("ggpubr")

## ----curateMGD, eval=FALSE----------------------------------------------------
#  library("curatedMetagenomicData")

## ----get_metadata, eval=FALSE-------------------------------------------------
#  meta.nielsen.full <- combined_metadata %>%
#      filter(dataset_name=='NielsenHB_2014')

## ----reason_meta_clean_2, eval=FALSE------------------------------------------
#  print(length(unique(meta.nielsen.full$subjectID)))
#  print(nrow(meta.nielsen.full))

## ----clean_metadata_1, eval=FALSE---------------------------------------------
#  meta.nielsen <- meta.nielsen.full %>%
#      select(sampleID, subjectID, study_condition, disease_subtype,
#          disease, age, country, number_reads, median_read_length, BMI) %>%
#      mutate(visit=str_extract(sampleID, '_[0-9]+$')) %>%
#      mutate(visit=str_remove(visit, '_')) %>%
#      mutate(visit=as.numeric(visit)) %>%
#      mutate(visit=case_when(is.na(visit)~0, TRUE~visit)) %>%
#      group_by(subjectID) %>%
#      filter(visit==min(visit)) %>%
#      ungroup() %>%
#      mutate(Sample_ID=sampleID) %>%
#      mutate(Group=case_when(disease=='healthy'~'CTR',
#                              TRUE~disease_subtype))

## ----clean_metadata_4, eval=FALSE---------------------------------------------
#  meta.nielsen <- meta.nielsen %>%
#      filter(Group %in% c('UC', 'CTR'))

## ----clean_metadata_3, eval=FALSE---------------------------------------------
#  meta.nielsen <- meta.nielsen %>%
#      mutate(Country=country)
#  meta.nielsen <- as.data.frame(meta.nielsen)
#  rownames(meta.nielsen) <- meta.nielsen$sampleID

## ----tax_profiles, eval=FALSE-------------------------------------------------
#  x <- 'NielsenHB_2014.metaphlan_bugs_list.stool'
#  feat <- curatedMetagenomicData(x=x, dryrun=FALSE)
#  feat <- feat[[x]]@assayData$exprs

## ----clean_tax_profiles, eval=FALSE-------------------------------------------
#  feat <- feat[grep(x=rownames(feat), pattern='s__'),]
#  feat <- feat[grep(x=rownames(feat),pattern='t__', invert = TRUE),]
#  feat <- t(t(feat)/100)

## ----clean_feature_names, eval=FALSE------------------------------------------
#  rownames(feat) <- str_extract(rownames(feat), 's__.*$')

## ----load_motus---------------------------------------------------------------
# base url for data download
data.loc <- 'https://zenodo.org/api/files/d81e429c-870f-44e0-a44a-2a4aa541b6c1/'
## metadata
meta.nielsen <- read_tsv(paste0(data.loc, 'meta_Nielsen.tsv'))
# also here, we have to remove repeated samplings and CD samples
meta.nielsen <- meta.nielsen %>%
    filter(Group %in% c('CTR', 'UC')) %>%
    group_by(Individual_ID) %>%
    filter(Sampling_day==min(Sampling_day)) %>%
    ungroup() %>%
    as.data.frame()
rownames(meta.nielsen) <- meta.nielsen$Sample_ID

## features
feat <- read.table(paste0(data.loc, 'metaHIT_motus.tsv'), 
                    stringsAsFactors = FALSE, sep='\t',
                    check.names = FALSE, quote = '', comment.char = '')
feat <- feat[,colSums(feat) > 0]
feat <- prop.table(as.matrix(feat), 2)

## ----siamcat------------------------------------------------------------------
# remove Danish samples
meta.nielsen.esp <- meta.nielsen[meta.nielsen$Country == 'ESP',]
sc.obj <- siamcat(feat=feat, meta=meta.nielsen.esp, label='Group', case='UC')

## ----feature_filtering--------------------------------------------------------
sc.obj <- filter.features(sc.obj, cutoff=1e-04,
                            filter.method = 'abundance')
sc.obj <- filter.features(sc.obj, cutoff=0.05,
                            filter.method='prevalence',
                            feature.type = 'filtered')

## ----assoc_plot, message=FALSE, warning=FALSE---------------------------------
sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha=0.1)
association.plot(sc.obj, fn.plot = './association_plot_nielsen.pdf', 
                panels = c('fc'))

## ----check_confounders, warning=FALSE-----------------------------------------
check.confounders(sc.obj, fn.plot = './confounders_nielsen.pdf')

## ----ml_workflow, eval=FALSE--------------------------------------------------
#  sc.obj <- normalize.features(sc.obj, norm.method = 'log.std',
#                              norm.param = list(log.n0=1e-06, sd.min.q=0))
#  ## Features normalized successfully.
#  sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5)
#  ## Features splitted for cross-validation successfully.
#  sc.obj <- train.model(sc.obj, method='lasso')
#  ## Trained lasso models successfully.
#  sc.obj <- make.predictions(sc.obj)
#  ## Made predictions successfully.
#  sc.obj <- evaluate.predictions(sc.obj)
#  ## Evaluated predictions successfully.

## ----model_eval_plot, eval=FALSE----------------------------------------------
#  model.evaluation.plot(sc.obj, fn.plot = './eval_plot_nielsen.pdf')
#  ## Plotted evaluation of predictions successfully to: ./eval_plot_nielsen.pdf

## ----model_interpretation_plot, eval=FALSE------------------------------------
#  model.interpretation.plot(sc.obj, consens.thres = 0.8,
#                              fn.plot = './interpret_nielsen.pdf')
#  ## Successfully plotted model interpretation plot to: ./interpret_nielsen.pdf

## ----confounders--------------------------------------------------------------
table(meta.nielsen$Group, meta.nielsen$Country)

## ----conf_start---------------------------------------------------------------
sc.obj.full <- siamcat(feat=feat, meta=meta.nielsen,
                        label='Group', case='UC')
sc.obj.full <- filter.features(sc.obj.full, cutoff=1e-04,
                                filter.method = 'abundance')
sc.obj.full <- filter.features(sc.obj.full, cutoff=0.05,
                                filter.method='prevalence',
                                feature.type = 'filtered')

## ----conf_country, eval=FALSE-------------------------------------------------
#  check.confounders(sc.obj.full, fn.plot = './confounders_dnk.pdf')

## ----assoc_plot_2, warning=FALSE, message=FALSE-------------------------------
sc.obj.full <- check.associations(sc.obj.full, log.n0 = 1e-06, alpha=0.1) 

## ----conf_assoc_plot, warning=FALSE-------------------------------------------
assoc.sp <- associations(sc.obj)
assoc.sp$species <- rownames(assoc.sp)
assoc.sp_dnk <- associations(sc.obj.full)
assoc.sp_dnk$species <- rownames(assoc.sp_dnk)

df.plot <- full_join(assoc.sp, assoc.sp_dnk, by='species')
df.plot %>%
    mutate(highlight=str_detect(species, 'formicigenerans')) %>%
    ggplot(aes(x=-log10(p.adj.x), y=-log10(p.adj.y), col=highlight)) +
    geom_point(alpha=0.3) +
        xlab('Spanish samples only\n-log10(q)') +
        ylab('Spanish and Danish samples only\n-log10(q)') +
        theme_classic() +
        theme(panel.grid.major = element_line(colour='lightgrey'),
            aspect.ratio = 1.3) +
        scale_colour_manual(values=c('darkgrey', '#D41645'), guide='none') +
        annotate('text', x=0.7, y=8, label='Dorea formicigenerans')

## ----dorea_plot---------------------------------------------------------------
# extract information out of the siamcat object
feat.dnk <- get.filt_feat.matrix(sc.obj.full)
label.dnk <- label(sc.obj.full)$label
country <- meta(sc.obj.full)$Country
names(country) <- rownames(meta(sc.obj.full))

df.plot <- tibble(dorea=log10(feat.dnk[
    str_detect(rownames(feat.dnk),'formicigenerans'),
    names(label.dnk)] + 1e-05),
    label=label.dnk, country=country) %>%
    mutate(label=case_when(label=='-1'~'CTR', TRUE~"UC")) %>%
    mutate(x_value=paste0(country, '_', label))

df.plot %>%
    ggplot(aes(x=x_value, y=dorea)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.08, stroke=0, alpha=0.2) +
        theme_classic() +
        xlab('') +
        ylab("log10(Dorea formicigenerans)") +
        stat_compare_means(comparisons = list(c('DNK_CTR', 'ESP_CTR'),
                                                c('DNK_CTR', 'ESP_UC'),
                                                c('ESP_CTR', 'ESP_UC')))

## ----ml_workflow_dnk, eval=FALSE----------------------------------------------
#  sc.obj.full <- normalize.features(sc.obj.full, norm.method = 'log.std',
#                                  norm.param = list(log.n0=1e-06, sd.min.q=0))
#  ## Features normalized successfully.
#  sc.obj.full <- create.data.split(sc.obj.full, num.folds = 5, num.resample = 5)
#  ## Features splitted for cross-validation successfully.
#  sc.obj.full <- train.model(sc.obj.full, method='lasso')
#  ## Trained lasso models successfully.
#  sc.obj.full <- make.predictions(sc.obj.full)
#  ## Made predictions successfully.
#  sc.obj.full <- evaluate.predictions(sc.obj.full)
#  ## Evaluated predictions successfully.

## ----eval_plot_comp, eval=FALSE-----------------------------------------------
#  model.evaluation.plot("Spanish samples only"=sc.obj,
#                      "Danish and Spanish samples"=sc.obj.full,
#                      fn.plot = './eval_plot_dnk.pdf')
#  ## Plotted evaluation of predictions successfully to: ./eval_plot_dnk.pdf

## ----ml_workflow_country, eval=FALSE------------------------------------------
#  meta.nielsen.country <- meta.nielsen[meta.nielsen$Group=='CTR',]
#  
#  sc.obj.country <- siamcat(feat=feat, meta=meta.nielsen.country,
#                              label='Country', case='ESP')
#  sc.obj.country <- filter.features(sc.obj.country, cutoff=1e-04,
#                              filter.method = 'abundance')
#  sc.obj.country <- filter.features(sc.obj.country, cutoff=0.05,
#                              filter.method='prevalence',
#                              feature.type = 'filtered')
#  sc.obj.country <- normalize.features(sc.obj.country, norm.method = 'log.std',
#                                      norm.param = list(log.n0=1e-06,
#                                          sd.min.q=0))
#  sc.obj.country <- create.data.split(sc.obj.country,
#                                      num.folds = 5, num.resample = 5)
#  sc.obj.country <- train.model(sc.obj.country, method='lasso')
#  sc.obj.country <- make.predictions(sc.obj.country)
#  sc.obj.country <- evaluate.predictions(sc.obj.country)
#  
#  print(eval_data(sc.obj.country)$auroc)
#  ## Area under the curve: 0.9701

## ----session_info-------------------------------------------------------------
sessionInfo()

