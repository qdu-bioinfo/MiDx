## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library("tidyverse")
library("SIAMCAT")

## ----load_data, message=FALSE-------------------------------------------------
# base url for data download
data.loc <- 'https://zenodo.org/api/files/d81e429c-870f-44e0-a44a-2a4aa541b6c1/'
# datasets
datasets <- c('metaHIT', 'Lewis_2015', 'He_2017', 'Franzosa_2019', 'HMP2')
# metadata
meta.all <- read_tsv(paste0(data.loc, 'meta_all_cd.tsv'))
# features
feat <- read.table(paste0(data.loc, 'feat_genus_cd.tsv'), 
                check.names = FALSE, stringsAsFactors = FALSE, quote = '', 
                sep='\t')
feat <- as.matrix(feat)
# check that metadata and features agree
stopifnot(all(colnames(feat) == meta.all$Sample_ID))

## ----table--------------------------------------------------------------------
table(meta.all$Study, meta.all$Group)

## ----meta_dereplicate---------------------------------------------------------
meta.ind <- meta.all %>% 
    group_by(Individual_ID) %>% 
    filter(Timepoint==min(Timepoint)) %>% 
    ungroup()

## ----compute_associations, message=FALSE, warning=FALSE-----------------------
assoc.list <- list()
for (d in datasets){
    # filter metadata and convert to dataframe
    meta.train <- meta.ind %>% 
        filter(Study==d) %>% 
        as.data.frame()
    rownames(meta.train) <- meta.train$Sample_ID

    # create SIAMCAT object
    sc.obj <- siamcat(feat=feat, meta=meta.train, label='Group', case='CD')
    # test for associations
    sc.obj <- check.associations(sc.obj, log.n0=1e-05, 
                                feature.type = 'original')
    # extract the associations and save them in the assoc.list
    temp <- associations(sc.obj)
    temp$genus <- rownames(temp)
    assoc.list[[d]] <- temp %>% 
        select(genus, fc, auc, p.adj) %>% 
        mutate(Study=d)
}
# combine all associations
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(genus!='unclassified')
head(df.assoc)

## ----genera_of_interest-------------------------------------------------------
genera.of.interest <- df.assoc %>% 
    group_by(genus) %>% 
    summarise(m=mean(auc), n.filt=any(auc < 0.25 | auc > 0.75), 
        .groups='keep') %>% 
    filter(n.filt) %>% 
    arrange(m)

## ----heatmap------------------------------------------------------------------
df.assoc %>% 
    # take only genera of interest
    filter(genus %in% genera.of.interest$genus) %>% 
    # convert to factor to enforce an ordering by mean AUC
    mutate(genus=factor(genus, levels = rev(genera.of.interest$genus))) %>% 
    # convert to factor to enforce ordering again
    mutate(Study=factor(Study, levels = datasets)) %>% 
    # annotate the cells in the heatmap with stars
    mutate(l=case_when(p.adj < 0.01~'*', TRUE~'')) %>%  
    ggplot(aes(y=genus, x=Study, fill=fc)) + 
        geom_tile() + 
        scale_fill_gradient2(low = '#3B6FB6', high='#D41645', mid = 'white', 
            limits=c(-2.7, 2.7), name='Generalized\nfold change') + 
        theme_minimal() + 
        geom_text(aes(label=l)) +
        theme(panel.grid = element_blank()) + 
        xlab('') + ylab('') +
        theme(axis.text = element_text(size=6))

## ----check_confounders, warning=FALSE-----------------------------------------
df.meta <- meta.ind %>% 
    as.data.frame()
rownames(df.meta) <- df.meta$Sample_ID
sc.obj <- siamcat(feat=feat, meta=df.meta, label='Group', case='CD')
check.confounders(sc.obj, fn.plot = './confounder_plot_cd_meta.pdf',
                feature.type='original')

## ----ml_meta_analysis, message=FALSE, warning=FALSE, eval=FALSE---------------
#  # create tibble to store all the predictions
#  auroc.all <- tibble(study.train=character(0),
#                      study.test=character(0),
#                      AUC=double(0))
#  # and a list to save the trained SIAMCAT objects
#  sc.list <- list()
#  for (i in datasets){
#      # restrict to a single study
#      meta.train <- meta.all %>%
#          filter(Study==i) %>%
#          as.data.frame()
#      rownames(meta.train) <- meta.train$Sample_ID
#  
#      ## take into account repeated sampling by including a parameters
#      ## in the create.data.split function
#      ## For studies with repeated samples, we want to block the
#      ## cross validation by the column 'Individual_ID'
#      block <- NULL
#      if (i %in% c('metaHIT', 'Lewis_2015', 'HMP2')){
#          block <- 'Individual_ID'
#          if (i == 'HMP2'){
#              # for the HMP2 dataset, the number of repeated sample per subject
#              # need to be reduced, because some subjects have been sampled
#              # 20 times, other only 5 times
#              meta.train <- meta.all %>%
#                  filter(Study=='HMP2') %>%
#                  group_by(Individual_ID) %>%
#                  sample_n(5, replace = TRUE) %>%
#                  distinct() %>%
#                  as.data.frame()
#              rownames(meta.train) <- meta.train$Sample_ID
#          }
#      }
#      # create SIAMCAT object
#      sc.obj.train <- siamcat(feat=feat, meta=meta.train,
#                              label='Group', case='CD')
#      # normalize features
#      sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
#          norm.param=list(log.n0=1e-05, sd.min.q=0),feature.type = 'original')
#      # Create data split
#      sc.obj.train <- create.data.split(sc.obj.train,
#          num.folds = 10, num.resample = 10, inseparable = block)
#      # train LASSO model
#      sc.obj.train <- train.model(sc.obj.train, method='lasso')
#  
#  
#      ## apply trained models to other datasets
#  
#      # loop through datasets again
#      for (i2 in datasets){
#          if (i == i2){
#              # make and evaluate cross-validation predictions (same dataset)
#              sc.obj.train <- make.predictions(sc.obj.train)
#              sc.obj.train <- evaluate.predictions(sc.obj.train)
#              auroc.all <- auroc.all %>%
#                  add_row(study.train=i, study.test=i,
#                      AUC=eval_data(sc.obj.train)$auroc %>% as.double())
#          } else {
#              # make and evaluate on the external datasets
#              # use meta.ind here, since we want only one sample per subject!
#              meta.test <- meta.ind %>%
#                  filter(Study==i2) %>%
#                  as.data.frame()
#              rownames(meta.test) <- meta.test$Sample_ID
#              sc.obj.test <- siamcat(feat=feat, meta=meta.test,
#                                      label='Group', case='CD')
#              # make holdout predictions
#              sc.obj.test <- make.predictions(sc.obj.train,
#                                              siamcat.holdout = sc.obj.test)
#              sc.obj.test <- evaluate.predictions(sc.obj.test)
#              auroc.all <- auroc.all %>%
#                  add_row(study.train=i, study.test=i2,
#                      AUC=eval_data(sc.obj.test)$auroc %>% as.double())
#          }
#      }
#      # save the trained model
#      sc.list[[i]] <- sc.obj.train
#  }

## ----load_aurocs, echo=FALSE, message=FALSE-----------------------------------
fn.in.auroc <- system.file(
    "extdata",
    "cd_meta_auroc.tsv",
    package = "SIAMCAT"
)
auroc.all <- read_tsv(fn.in.auroc)

## ----get_test_average---------------------------------------------------------
test.average <- auroc.all %>% 
    filter(study.train!=study.test) %>% 
    group_by(study.test) %>% 
    summarise(AUC=mean(AUC), .groups='drop') %>% 
    mutate(study.train="Average")

## ----plot_auroc---------------------------------------------------------------
# combine AUROC values with test average
bind_rows(auroc.all, test.average) %>% 
    # highlight cross validation versus transfer results
    mutate(CV=study.train == study.test) %>%
    # for facetting later
    mutate(split=case_when(study.train=='Average'~'Average', TRUE~'none')) %>% 
    mutate(split=factor(split, levels = c('none', 'Average'))) %>% 
    # convert to factor to enforce ordering
    mutate(study.train=factor(study.train, levels=c(datasets, 'Average'))) %>% 
    mutate(study.test=factor(study.test, 
                            levels=c(rev(datasets),'Average'))) %>% 
    ggplot(aes(y=study.test, x=study.train, fill=AUC, size=CV, color=CV)) +
        geom_tile() + theme_minimal() +
        # text in tiles
        geom_text(aes_string(label="format(AUC, digits=2)"), 
            col='white', size=2)+
        # color scheme
        scale_fill_gradientn(colours=rev(c('darkgreen','forestgreen', 
                                        'chartreuse3','lawngreen', 
                                        'yellow')), limits=c(0.5, 1)) +
        # axis position/remove boxes/ticks/facet background/etc.
        scale_x_discrete(position='top') + 
        theme(axis.line=element_blank(), 
                axis.ticks = element_blank(), 
                axis.text.x.top = element_text(angle=45, hjust=.1), 
                panel.grid=element_blank(), 
                panel.border=element_blank(), 
                strip.background = element_blank(), 
                strip.text = element_blank()) + 
        xlab('Training Set') + ylab('Test Set') + 
        scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
        scale_size_manual(values=c(0, 1), guide=FALSE) + 
        facet_grid(~split, scales = 'free', space = 'free')

## ----weights, warning=FALSE, message=FALSE, eval=FALSE------------------------
#  weight.list <- list()
#  for (d in datasets){
#      sc.obj.train <- sc.list[[d]]
#      # extract the feature weights out of the SIAMCAT object
#      temp <- feature_weights(sc.obj.train)
#      temp$genus <- rownames(temp)
#      # save selected info in the weight.list
#      weight.list[[d]] <- temp %>%
#          select(genus, median.rel.weight, mean.rel.weight, percentage) %>%
#          mutate(Study=d) %>%
#          mutate(r.med=rank(-abs(median.rel.weight)),
#              r.mean=rank(-abs(mean.rel.weight)))
#  }
#  # combine all feature weights into a single tibble
#  df.weights <- bind_rows(weight.list)
#  df.weights <- df.weights %>% filter(genus!='unclassified')

## ----load_weights, echo=FALSE, message=FALSE----------------------------------
fn.in.weights<- system.file(
    "extdata",
    "cd_meta_weights.tsv",
    package = "SIAMCAT"
)
df.weights <- read_tsv(fn.in.weights)

## ----plot_weights_heatmap, warning=FALSE--------------------------------------
# compute absolute feature weights
abs.weights <- df.weights %>% 
    group_by(Study) %>% 
    summarise(sum.median=sum(abs(median.rel.weight)),
                sum.mean=sum(abs(mean.rel.weight)),
                .groups='drop')

df.weights %>% 
    full_join(abs.weights) %>% 
    # normalize by the absolute model size
    mutate(median.rel.weight=median.rel.weight/sum.median) %>% 
    # only include genera of interest
    filter(genus %in% genera.of.interest$genus) %>% 
    # highlight feature rank for the top 20 features
    mutate(r.med=case_when(r.med > 20~NA_real_, TRUE~r.med)) %>%
    # enforce the correct ordering by converting to factors again
    mutate(genus=factor(genus, levels = rev(genera.of.interest$genus))) %>% 
    mutate(Study=factor(Study, levels = datasets)) %>% 
    ggplot(aes(y=genus, x=Study, fill=median.rel.weight)) + 
        geom_tile() + 
        scale_fill_gradientn(colours=rev(
            c('#007A53', '#009F4D', "#6CC24A", 'white',
            "#EFC06E", "#FFA300", '#BE5400')), 
            limits=c(-0.15, 0.15)) +
        theme_minimal() + 
        geom_text(aes(label=r.med), col='black', size= 2) +
        theme(panel.grid = element_blank()) + 
        xlab('') + ylab('') +
        theme(axis.text = element_text(size=6))



## ----session_info-------------------------------------------------------------
sessionInfo()

