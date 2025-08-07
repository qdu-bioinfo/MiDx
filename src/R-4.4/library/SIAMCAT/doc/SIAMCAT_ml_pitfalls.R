## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE)

## ----load, message=FALSE, warning=FALSE---------------------------------------
library("tidyverse")
library("SIAMCAT")

## ----load_curatedMetagenomicsData, eval=FALSE---------------------------------
#  library("curatedMetagenomicData")

## ----load_data_thomas, eval=FALSE---------------------------------------------
#  x <- 'ThomasAM_2018a.metaphlan_bugs_list.stool'
#  feat.t <- curatedMetagenomicData(x=x, dryrun=FALSE)
#  feat.t <- feat.t[[x]]@assayData$exprs
#  # clean up metaphlan profiles to contain only species-level abundances
#  feat.t <- feat.t[grep(x=rownames(feat.t), pattern='s__'),]
#  feat.t <- feat.t[grep(x=rownames(feat.t),pattern='t__', invert = TRUE),]
#  stopifnot(all(colSums(feat.t) != 0))
#  feat.t <- t(t(feat.t)/100)

## ----load_data_zeller, eval=FALSE---------------------------------------------
#  x <- 'ZellerG_2014.metaphlan_bugs_list.stool'
#  feat.z <- curatedMetagenomicData(x=x, dryrun=FALSE)
#  feat.z <- feat.z[[x]]@assayData$exprs
#  # clean up metaphlan profiles to contain only species-level abundances
#  feat.z <- feat.z[grep(x=rownames(feat.z), pattern='s__'),]
#  feat.z <- feat.z[grep(x=rownames(feat.z),pattern='t__', invert = TRUE),]
#  stopifnot(all(colSums(feat.z) != 0))
#  feat.z <- t(t(feat.z)/100)

## ----metadata, eval=FALSE-----------------------------------------------------
#  meta.t <- combined_metadata %>%
#      filter(dataset_name == 'ThomasAM_2018a') %>%
#      filter(study_condition %in% c('control', 'CRC'))
#  rownames(meta.t) <- meta.t$sampleID
#  meta.z <- combined_metadata %>%
#      filter(dataset_name == 'ZellerG_2014') %>%
#      filter(study_condition %in% c('control', 'CRC'))
#  rownames(meta.z) <- meta.z$sampleID

## ----combine features, eval=FALSE---------------------------------------------
#  species.union <- union(rownames(feat.t), rownames(feat.z))
#  # add Zeller_2014-only species to the Thomas_2018 matrix
#  add.species <- setdiff(species.union, rownames(feat.t))
#  feat.t <- rbind(feat.t,
#              matrix(0, nrow=length(add.species), ncol=ncol(feat.t),
#                  dimnames = list(add.species, colnames(feat.t))))
#  
#  # add Thomas_2018-only species to the Zeller_2014 matrix
#  add.species <- setdiff(species.union, rownames(feat.z))
#  feat.z <- rbind(feat.z,
#              matrix(0, nrow=length(add.species), ncol=ncol(feat.z),
#                  dimnames = list(add.species, colnames(feat.z))))

## ----setup_2------------------------------------------------------------------
fs.cutoff <- c(20, 100, 250)

auroc.all <- tibble(cutoff=character(0), type=character(0), 
                    study.test=character(0), AUC=double(0))

## ----train_full, eval=FALSE---------------------------------------------------
#  sc.obj.t <- siamcat(feat=feat.t, meta=meta.t,
#                      label='study_condition', case='CRC')
#  sc.obj.t <- filter.features(sc.obj.t, filter.method = 'prevalence',
#                              cutoff = 0.01)
#  sc.obj.t <- normalize.features(sc.obj.t, norm.method = 'log.std',
#                                  norm.param=list(log.n0=1e-05, sd.min.q=0))
#  sc.obj.t <- create.data.split(sc.obj.t,
#                                  num.folds = 10, num.resample = 10)
#  sc.obj.t <- train.model(sc.obj.t, method='lasso')
#  sc.obj.t <- make.predictions(sc.obj.t)
#  sc.obj.t <- evaluate.predictions(sc.obj.t)
#  
#  auroc.all <- auroc.all %>%
#      add_row(cutoff='full', type='correct',
#              study.test='Thomas_2018',
#              AUC=as.numeric(sc.obj.t@eval_data$auroc)) %>%
#      add_row(cutoff='full', type='incorrect', study.test='Thomas_2018',
#              AUC=as.numeric(sc.obj.t@eval_data$auroc))

## ----ext_val_full, eval=FALSE-------------------------------------------------
#  sc.obj.z <- siamcat(feat=feat.z, meta=meta.z,
#                      label='study_condition', case='CRC')
#  sc.obj.z <- make.predictions(sc.obj.t, sc.obj.z)
#  sc.obj.z <- evaluate.predictions(sc.obj.z)
#  auroc.all <- auroc.all %>%
#      add_row(cutoff='full', type='correct',
#              study.test='Zeller_2014',
#              AUC=as.numeric(sc.obj.z@eval_data$auroc)) %>%
#      add_row(cutoff='full', type='incorrect',
#              study.test='Zeller_2014',
#              AUC=as.numeric(sc.obj.z@eval_data$auroc))

## ----train_global, eval=FALSE-------------------------------------------------
#  sc.obj.t <- check.associations(sc.obj.t, detect.lim = 1e-05,
#                                  fn.plot = 'assoc_plot.pdf')
#  mat.assoc <- associations(sc.obj.t)
#  mat.assoc$species <- rownames(mat.assoc)
#  # sort by p-value
#  mat.assoc <- mat.assoc %>% as_tibble() %>% arrange(p.val)

## ----train_global_2, eval=FALSE-----------------------------------------------
#  for (x in fs.cutoff){
#      # select x number of features based on p-value ranking
#      feat.train.red <- feat.t[mat.assoc %>%
#                                  slice(seq_len(x)) %>%
#                                  pull(species),]
#      sc.obj.t.fs <- siamcat(feat=feat.train.red, meta=meta.t,
#                              label='study_condition', case='CRC')
#      # normalize the features without filtering
#      sc.obj.t.fs <- normalize.features(sc.obj.t.fs, norm.method = 'log.std',
#          norm.param=list(log.n0=1e-05,sd.min.q=0),feature.type = 'original')
#      # take the same cross validation split as before
#      data_split(sc.obj.t.fs) <- data_split(sc.obj.t)
#      # train
#      sc.obj.t.fs <- train.model(sc.obj.t.fs, method = 'lasso')
#      # make predictions
#      sc.obj.t.fs <- make.predictions(sc.obj.t.fs)
#      # evaluate predictions and record the result
#      sc.obj.t.fs <- evaluate.predictions(sc.obj.t.fs)
#      auroc.all <- auroc.all %>%
#          add_row(cutoff=as.character(x), type='incorrect',
#                  study.test='Thomas_2018',
#                  AUC=as.numeric(sc.obj.t.fs@eval_data$auroc))
#      # apply to the external dataset and record the result
#      sc.obj.z <- siamcat(feat=feat.z, meta=meta.z,
#                          label='study_condition', case='CRC')
#      sc.obj.z <- make.predictions(sc.obj.t.fs, sc.obj.z)
#      sc.obj.z <- evaluate.predictions(sc.obj.z)
#      auroc.all <- auroc.all %>%
#          add_row(cutoff=as.character(x), type='incorrect',
#                  study.test='Zeller_2014',
#                  AUC=as.numeric(sc.obj.z@eval_data$auroc))
#  }

## ----train_nested, eval=FALSE-------------------------------------------------
#  for (x in fs.cutoff){
#      # train using the original SIAMCAT object
#      # with correct version of feature selection
#      sc.obj.t.fs <- train.model(sc.obj.t, method = 'lasso', perform.fs = TRUE,
#          param.fs = list(thres.fs = x,method.fs = "AUC",direction='absolute'))
#      # make predictions
#      sc.obj.t.fs <- make.predictions(sc.obj.t.fs)
#      # evaluate predictions and record the result
#      sc.obj.t.fs <- evaluate.predictions(sc.obj.t.fs)
#      auroc.all <- auroc.all %>%
#          add_row(cutoff=as.character(x), type='correct',
#                  study.test='Thomas_2018',
#                  AUC=as.numeric(sc.obj.t.fs@eval_data$auroc))
#      # apply to the external dataset and record the result
#      sc.obj.z <- siamcat(feat=feat.z, meta=meta.z,
#                          label='study_condition', case='CRC')
#      sc.obj.z <- make.predictions(sc.obj.t.fs, sc.obj.z)
#      sc.obj.z <- evaluate.predictions(sc.obj.z)
#      auroc.all <- auroc.all %>%
#          add_row(cutoff=as.character(x), type='correct',
#                  study.test='Zeller_2014',
#                  AUC=as.numeric(sc.obj.z@eval_data$auroc))
#  }

## ----data, echo=FALSE---------------------------------------------------------
auroc.all <- tibble(
    cutoff=rep(rep(c('20', '100', '250', 'full'), each=2), 2),
    type=rep(c('incorrect', 'correct'), 8),
    study.test=rep(c('Thomas_2018', 'Zeller_2014'), each=8),
    AUC=c(0.809, 0.608, 0.812, 0.659, 0.727, 0.678, 0.677, 0.677,
            0.620, 0.688, 0.694, 0.732, 0.737, 0.737, 0.736, 0.736))

## ----plot_auroc---------------------------------------------------------------
auroc.all %>%
    # facetting for plotting
    mutate(split=case_when(study.test=="Thomas_2018"~
                            'Cross validation (Thomas 2018)',
                        TRUE~"External validation (Zeller 2014)")) %>%
    # convert to factor to enforce ordering
    mutate(cutoff=factor(cutoff, levels = c(fs.cutoff, 'full'))) %>%
    ggplot(aes(x=cutoff, y=AUC, col=type)) +
        geom_point() + geom_line(aes(group=type)) +
        facet_grid(~split) +
        scale_y_continuous(limits = c(0.5, 1), expand = c(0,0)) +
        xlab('Features selected') +
        ylab('AUROC') +
        theme_bw() + 
        scale_colour_manual(values = c('correct'='blue', 'incorrect'='red'),
            name='Feature selection procedure') + 
        theme(panel.grid.minor = element_blank(), legend.position = 'bottom')
    

## ----load_data_ibd------------------------------------------------------------
data.loc <- 'https://zenodo.org/api/files/d81e429c-870f-44e0-a44a-2a4aa541b6c1/'

# metadata
meta.all <- read_tsv(paste0(data.loc, 'meta_all_cd.tsv'))

# features
feat.motus <- read.table(paste0(data.loc, 'feat_rel_filt_cd.tsv'),
                        sep='\t', stringsAsFactors = FALSE,
                        check.names = FALSE)

## ----no_samples_per_indiv-----------------------------------------------------
x <- meta.all %>% 
    group_by(Study, Group) %>% 
    summarise(n.all=n(), .groups='drop')
y <- meta.all %>% 
    select(Study, Group, Individual_ID) %>% 
    distinct() %>% 
    group_by(Study, Group) %>% 
    summarize(n.indi=n(),  .groups='drop')
full_join(x,y)

## ----hmp_samples--------------------------------------------------------------
meta.all %>% 
    filter(Study=='HMP2') %>% 
    group_by(Individual_ID) %>% 
    summarise(n=n(), .groups='drop') %>% 
    pull(n) %>% hist(20)

# sample 5 samples per individual
meta.train <- meta.all %>% 
    filter(Study=='HMP2') %>% 
    group_by(Individual_ID) %>%
    sample_n(5, replace = TRUE) %>%
    distinct() %>%
    as.data.frame()
rownames(meta.train) <- meta.train$Sample_ID

## ----meta_ind-----------------------------------------------------------------
meta.ind <- meta.all %>% 
    group_by(Individual_ID) %>% 
    filter(Timepoint==min(Timepoint)) %>% 
    ungroup()

## ----create_auc_tibble--------------------------------------------------------
auroc.all <- tibble(type=character(0), study.test=character(0), AUC=double(0))

## ----train_incorrect_ibd_models, eval=FALSE-----------------------------------
#  sc.obj <- siamcat(feat=feat.motus, meta=meta.train,
#                      label='Group', case='CD')
#  sc.obj <- normalize.features(sc.obj, norm.method = 'log.std',
#      norm.param=list(log.n0=1e-05,sd.min.q=1),feature.type = 'original')
#  sc.obj.naive <- create.data.split(sc.obj, num.folds = 10, num.resample = 10)
#  sc.obj.naive <- train.model(sc.obj.naive, method='lasso')
#  sc.obj.naive <- make.predictions(sc.obj.naive)
#  sc.obj.naive <- evaluate.predictions(sc.obj.naive)
#  auroc.all <- auroc.all %>%
#      add_row(type='naive', study.test='HMP2',
#          AUC=as.numeric(eval_data(sc.obj.naive)$auroc))

## ----train_correct_ibd_models, eval=FALSE-------------------------------------
#  sc.obj.block <- create.data.split(sc.obj, num.folds = 10, num.resample = 10,
#                                  inseparable = 'Individual_ID')
#  sc.obj.block <- train.model(sc.obj.block, method='lasso')
#  sc.obj.block <- make.predictions(sc.obj.block)
#  sc.obj.block <- evaluate.predictions(sc.obj.block)
#  auroc.all <- auroc.all %>%
#      add_row(type='blocked', study.test='HMP2',
#          AUC=as.numeric(eval_data(sc.obj.block)$auroc))

## ----apply_ibd_models, eval=FALSE---------------------------------------------
#  for (i in setdiff(unique(meta.all$Study), 'HMP2')){
#      meta.test <- meta.ind %>%
#          filter(Study==i) %>%
#          as.data.frame()
#      rownames(meta.test) <- meta.test$Sample_ID
#      # apply naive model
#      sc.obj.test <- siamcat(feat=feat.motus, meta=meta.test,
#                              label='Group', case='CD')
#      sc.obj.test <- make.predictions(sc.obj.naive, sc.obj.test)
#      sc.obj.test <- evaluate.predictions(sc.obj.test)
#      auroc.all <- auroc.all %>%
#      add_row(type='naive', study.test=i,
#              AUC=as.numeric(eval_data(sc.obj.test)$auroc))
#      # apply blocked model
#      sc.obj.test <- siamcat(feat=feat.motus, meta=meta.test,
#                              label='Group', case='CD')
#      sc.obj.test <- make.predictions(sc.obj.block, sc.obj.test)
#      sc.obj.test <- evaluate.predictions(sc.obj.test)
#      auroc.all <- auroc.all %>%
#          add_row(type='blocked', study.test=i,
#                  AUC=as.numeric(eval_data(sc.obj.test)$auroc))
#  }

## ----load_results_dp, echo=FALSE----------------------------------------------
auroc.all <- tibble(
    type=rep(c('naive', 'blocked'), 5),
    study.test=(rep(c('metaHIT', 'Lewis_2015', 'He_2017',
        'Franzosa_2019', 'HMP2'), each=2)),
    AUC=c(0.77, 0.82, 0.80, 0.82, 0.788, 0.855, 0.739, 0.774, 0.988, 0.667))

## ----plot_results_db----------------------------------------------------------
auroc.all %>%
    # convert to factor to enforce ordering
    mutate(type=factor(type, levels = c('naive', 'blocked'))) %>%
    # facetting for plotting
    mutate(CV=case_when(study.test=='HMP2'~'CV', 
                        TRUE~'External validation')) %>%
    ggplot(aes(x=study.test, y=AUC, fill=type)) +
        geom_bar(stat='identity', position = position_dodge(), col='black') +
        theme_bw() +
        coord_cartesian(ylim=c(0.5, 1)) +
        scale_fill_manual(values=c('red', 'blue'), name='') +
        facet_grid(~CV, space = 'free', scales = 'free') +
        xlab('') + ylab('AUROC') +
        theme(legend.position = c(0.8, 0.8))

## ----session_info-------------------------------------------------------------
sessionInfo()

