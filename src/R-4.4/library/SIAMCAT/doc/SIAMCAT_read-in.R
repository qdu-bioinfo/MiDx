## ----feat_file, message=FALSE-------------------------------------------------
library(SIAMCAT)
fn.in.feat  <- system.file(
    "extdata",
    "feat_crc_zeller_msb_mocat_specI.tsv",
    package = "SIAMCAT"
)

## ----read_feat, message=FALSE-------------------------------------------------
feat <- read.table(fn.in.feat, sep='\t',
    header=TRUE, quote='',
    stringsAsFactors = FALSE, check.names = FALSE)
# look at some features
feat[110:114, 1:2]

## ----meta_file, message=FALSE-------------------------------------------------
fn.in.meta  <- system.file(
    "extdata",
    "num_metadata_crc_zeller_msb_mocat_specI.tsv",
    package = "SIAMCAT"
)

## ----read_meta, warning=FALSE-------------------------------------------------
meta <- read.table(fn.in.meta, sep='\t',
    header=TRUE, quote='',
    stringsAsFactors = FALSE, check.names = FALSE)
head(meta)

## ----create_label, results="hide", warning=FALSE, eval=FALSE------------------
#  label <- create.label(meta=meta, label="diagnosis",
#      case = 1, control=0)

## ----create_label_2, warning=FALSE--------------------------------------------
label <- create.label(meta=meta, label="diagnosis",
    case = 1, control=0,
    p.lab = 'cancer', n.lab = 'healthy')
label$info

## ----lefse_file, message=FALSE------------------------------------------------
fn.in.lefse<- system.file(
    "extdata",
    "LEfSe_crc_zeller_msb_mocat_specI.tsv",
    package = "SIAMCAT"
)

## ----read_lefse_file, results="hide", warning=FALSE---------------------------
meta.and.features <- read.lefse(fn.in.lefse,
    rows.meta = 1:6, row.samples = 7)
meta <- meta.and.features$meta
feat <- meta.and.features$feat

## ----lefse_label, results="hide", warning=FALSE-------------------------------
label <- create.label(meta=meta, label="label", case = "cancer")

## ----read_metagenome_seq, results="hide", warning=FALSE, eval=FALSE-----------
#  fn.in.feat  <- system.file(
#      "extdata",
#      "CHK_NAME.otus.count.csv",
#      package = "metagenomeSeq"
#  )
#  feat <- read.table(fn.in.feat, sep='\t',
#      header=TRUE, quote='', row.names = 1,
#      stringsAsFactors = FALSE, check.names = FALSE
#  )

## ----create_from_phyloseq, results="hide", warning=FALSE, eval=TRUE-----------
data("GlobalPatterns") ## phyloseq example data
label <- create.label(meta=sample_data(GlobalPatterns),
    label = "SampleType",
    case = c("Freshwater", "Freshwater (creek)", "Ocean"))
# run the constructor function
siamcat <- siamcat(phyloseq=GlobalPatterns, label=label)

## ----constructor, eval=FALSE--------------------------------------------------
#  siamcat <- siamcat(feat=feat, label=label, meta=meta)

## ----constructor_phyloseq, eval=FALSE-----------------------------------------
#  siamcat <- siamcat(phyloseq=phyloseq, label=label)

## ----eval=FALSE---------------------------------------------------------------
#  eval_data(siamcat)

## ----eval=FALSE---------------------------------------------------------------
#  label(siamcat) <- new_label

## -----------------------------------------------------------------------------
phyloseq <- physeq(siamcat)
tax_tab <- tax_table(phyloseq)
head(tax_tab)

## -----------------------------------------------------------------------------
sessionInfo()

