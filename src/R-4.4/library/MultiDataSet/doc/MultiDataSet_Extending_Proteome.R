## ----create_data, eval = TRUE, echo = FALSE-----------------------------------
#setwd("/home/kuragari/Downloads/use case")
#setwd("C:/Users/chernandez/shared/Projects/multidataset/Tutorials")
fdata <- data.frame(
    feature = c("Adiponectin", "CRP", "APO.A1", "APO.B", "APO.E"),
    LoD.T = c(4.0000E+07, 6.8333E+06, 1.6250E+09, 7.5000E+08, 3.0556E+08),
    LoD.B = c(1.6461E+05, 2.8121E+04,	2.0062E+07,	9.2593E+06,	3.7723E+06),
    unit = c("pg/mL", "pg/mL", "pg/mL", "pg/mL", "pg/mL")
)


pc <- textConnection("sample,gender,plate,kit
sp001,male,P01,5
sp002,male,P01,5
sp003,female,P01,5
sp004,male,P01,5
sp005,female,P01,5
sp006,female,P01,5
sp007,male,P01,5
sp009,male,P01,5
sp010,male,P01,5
sp012,male,P01,5
sp013,male,P01,5
sp014,female,P01,5
sp015,male,P01,5
sp016,male,P01,5
sp017,male,P01,5
sp018,male,P01,5
sp019,male,P01,5
sp020,male,P01,5
sp021,female,P01,5
sp022,male,P01,5
sp023,male,P01,5
sp024,female,P01,5
sp025,male,P01,5
sp026,male,P01,5
sp027,female,P01,5
sp028,male,P01,5
sp029,female,P01,5
sp030,female,P01,5
sp031,male,P01,5
sp032,female,P01,5
sp033,male,P01,5
sp034,male,P01,5
sp035,female,P01,5")

pdata <- read.csv(pc)

ac <- textConnection("sample,Adiponectin,CRP,APO.A1,APO.B,APO.E
sp001,15078717.84,256079.1978,166102189.6,81989413.94,12088576.71
sp002,13216433.82,59688.53882,122910159.9,75276925.92,8662032.005
sp003,25639817.37,272556.6707,132103162.9,58156884.1,11317489.32
sp004,28042622.37,981180.998,163470320.6,72249107.96,9534858.626
sp005,37507509.6,1434113.656,186415820,85977072.27,20772027.96
sp006,27201037.83,7140776.893,180018819.1,115016660.7,19020968.34
sp007,19734082.42,849104.2897,95178150.98,74531939.7,12575510.74
sp009,24337465.91,2083163.005,119909512.6,75162215.53,13472623.87
sp010,29346604.23,387120.4506,122159943.6,79429725.86,16315028.32
sp012,23471847.37,313730.4201,162530447.6,61607843.14,6950787.103
sp013,26464829.9,231195.7943,182276383.6,91789932.63,10467398.75
sp014,15188363.38,13795018.08,179266346.9,142227471.5,10822788.47
sp015,28585011.4,48724.14524,182652665.8,98888608.08,7252478.096
sp016,19105671.52,110506.2048,140362398.3,64931137.97,8072973.965
sp017,26570792.93,991854.9577,145432469.7,80533037.64,16914321.56
sp018,23968478.51,234826.4158,127787432.4,65209657.88,7924737.258
sp019,22881995.44,9321008.038,144681259.3,111528556.5,8072973.965
sp020,23118772.5,16891530.05,196202782.9,80823891.6,14734492.79
sp021,22942701.12,139509.27,92932170.07,59302740.96,4131144.701
sp022,25317279.56,103616.9905,115034724.1,75678689.62,7252478.096
sp023,26502576.75,79198.05423,209006876.8,110535951,20643082.4
sp024,25540633.93,334148.4346,182276383.6,84506590.67,14632785.98
sp025,34258505.11,654613.9767,115034724.1,90954853.78,18563191.06
sp026,17832612.2,229986.4403,111660816.1,61993740.5,10964463.6
sp027,20439508.67,534718.8069,157455877.1,56311320.18,7327611.264
sp028,31333461.27,179608.5461,113535112.2,69359486.56,12852515.92
sp029,18419548.11,271333.7529,100420374.9,41698956.49,10110231.2
sp030,15699158.6,1011889.525,147310628.7,54048254.2,8735244.722
sp031,20726453.63,479871.0422,98922372.45,58156884.1,10503016.43
sp032,33833352.39,532158.9756,214657722.7,84800277.57,19866893.4
sp033,14598125.04,179608.5461,84887644.72,86212822.7,13952043.6
sp034,24457846.67,1219055.651,124410700.1,85918155,18235193.56
sp035,11560113.29,134838.8042,92557881.06,34647582.51,6950787.103")

adata <- read.csv(ac)
rm(pc, ac)


write.table(fdata, file="feature_data.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(pdata, file="pheno_data.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(adata, file="assay_data.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
rm(fdata, pdata, adata)

## ----Define Proteome, message=FALSE-------------------------------------------
library(Biobase)
setClass (
    Class = "ProteomeSet",
    contains = "eSet", 
    prototype = prototype(new("VersionedBiobase",
                              versions = c(classVersion("eSet"), ProteomeSet = "1.0.0")))
)

## ----Define Proteome Validity, message=FALSE, results="hide"------------------
setValidity("ProteomeSet", function(object) {
    
    ## Check that object has the slots 'prots' and 'raw' in assayData
    msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("prots", "raw")))
    
    ## Check that objects has the columns 'LoD.T' and 'LoD.B' in featureData
    msg <- validMsg(msg, ifelse(any(!c("LoD.T", "LoD.B") %in% fvarLabels(object)), "fData must contain columns 'LoD.T' and 'LoD.B'", FALSE))
    if (is.null(msg)){
        TRUE
    } else{
        msg
    }
}) 

## ----read_ldset---------------------------------------------------------------
read_ldset <- function(assayFile, phenoFile, featureFile, sep="\t") {
    ## Load the threes files that will be used to create the ProteomeSet
    adata <- read.delim(assayFile, sep=sep, header=TRUE, row.names="sample")
    pdata <- read.delim(phenoFile, sep=sep, header=TRUE, row.names="sample")
    fdata <- read.delim(featureFile, sep=sep, header=TRUE, row.names="feature")
    ## /
    
    ## Check that proteins in assay data are the same in feature data
    if(!identical(colnames(adata), rownames(fdata))) {
        stop("Features names in assay data (columns) are not equal to ",
             "features names in feature data (rownames).")
    }
    ##/
    
    ## Check that feature data include columns LoD.B and LoD.T
    if(sum(c("LoD.T", "LoD.B") %in% colnames(fdata)) != 2) {
        stop("Feature data must contain two columns labeled 'LoD.T' (top ",
             "limit of dectection) and 'LoD.B (bottom limit of dectection)")
    }
    ## /
    
    ## Perform the transformation of the protein level of expression
    low <- fdata[colnames(adata), "LoD.B"]
    up <- fdata[colnames(adata), "LoD.T"]
    names(low) <- names(up) <- rownames(fdata)
    faux <- function(x, low, up) {
        x[x < low] <- as.double(low / 2)
        x[x > up] <- as.double(up * 1.5)
        x
    }
    tadata <- mapply(FUN = faux, x = as.list(adata), low = as.list(low), up = as.list(up))
    dimnames(tadata) <- dimnames(adata)
    ## /
    
    ## Create the ExpressionSet with the two matrices
    prot <- new("ProteomeSet",
                assayData = assayDataNew("environment", prots = t(tadata), raw = t(adata)),
                phenoData = AnnotatedDataFrame(pdata),
                featureData = AnnotatedDataFrame(fdata)
    )
    ## /
    
    ## Check that the new ProteomeSet is valid
    validObject(prot)
    ## /
    
    return(prot)
}

## ----create_generic-----------------------------------------------------------
setGeneric("add_prot", function(object, protSet, warnings = TRUE, ...)
    standardGeneric("add_prot")
)

## ----message = FALSE----------------------------------------------------------
library(MultiDataSet)
setMethod(
    f = "add_prot",
    signature = c("MultiDataSet", "ProteomeSet"),
    definition = function(object, protSet, warnings = TRUE, ...) {
        ## Add given ProteomeSet as 'proteome'
        object <- MultiDataSet::add_eset(object, protSet, dataset.type = "proteome", GRanges = NA, ...)
        ## /
        return(object)
    }
)

## ----load_ps------------------------------------------------------------------
## Create a ProteomeSet with protein data
ps <- read_ldset(assayFile="assay_data.tsv", 
                 phenoFile="pheno_data.tsv",
                 featureFile="feature_data.tsv"
)
ps

## ----add_prot_md_1, warning=FALSE---------------------------------------------
md <- createMultiDataSet()
md <- add_prot(md, ps)

## ----show_md_1----------------------------------------------------------------
names(md)

## ----show_md_2----------------------------------------------------------------
md

