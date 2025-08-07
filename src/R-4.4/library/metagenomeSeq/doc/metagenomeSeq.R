## ----include=FALSE---------------------------------------------
require(knitr)
opts_chunk$set(concordance=TRUE,tidy=TRUE)

## ----config,echo=FALSE------------------------------------
options(width = 60)
options(continue=" ")
options(warn=-1)
set.seed(42)

## ----requireMetagenomeSeq,warning=FALSE,message=FALSE-----
library(metagenomeSeq)

## ----loadBiom---------------------------------------------
# reading in a biom file
library(biomformat)
biom_file <- system.file("extdata", "min_sparse_otu_table.biom", package = "biomformat")
b <- read_biom(biom_file)
biom2MRexperiment(b)

## ----writeBiom,eval=FALSE---------------------------------
# data(mouseData)
# # options include to normalize or not
# b <- MRexperiment2biom(mouseData)
# write_biom(b,biom_file="~/Desktop/otu_table.biom")

## ----loadData---------------------------------------------
dataDirectory <- system.file("extdata", package="metagenomeSeq")
lung = loadMeta(file.path(dataDirectory,"CHK_NAME.otus.count.csv")) 
dim(lung$counts)

## ----loadTaxa---------------------------------------------
taxa = read.delim(file.path(dataDirectory,"CHK_otus.taxonomy.csv"),stringsAsFactors=FALSE)

## ----loadClin---------------------------------------------
clin = loadPhenoData(file.path(dataDirectory,"CHK_clinical.csv"),tran=TRUE)
ord = match(colnames(lung$counts),rownames(clin)) 
clin = clin[ord,]
head(clin[1:2,])

## ----createMRexperiment1----------------------------------
phenotypeData = AnnotatedDataFrame(clin)
phenotypeData

## ----createMRexperiment2----------------------------------
OTUdata = AnnotatedDataFrame(taxa)
OTUdata

## ----createMRexperiment3,tidy=FALSE-----------------------
obj = newMRexperiment(lung$counts,phenoData=phenotypeData,featureData=OTUdata)
# Links to a paper providing further details can be included optionally.
# experimentData(obj) = annotate::pmid2MIAME("21680950")
obj

## ----dataset1,tidy=FALSE----------------------------------
data(lungData)
lungData

## ----dataset2,tidy=FALSE----------------------------------
data(mouseData)
mouseData

## ----pdata------------------------------------------------
phenoData(obj)
head(pData(obj),3)

## ----fdata------------------------------------------------
featureData(obj)
head(fData(obj)[,-c(2,10)],3)

## ----MRcounts---------------------------------------------
head(MRcounts(obj[,1:2]))

## ---------------------------------------------------------
featuresToKeep = which(rowSums(obj)>=100)
samplesToKeep = which(pData(obj)$SmokingStatus=="Smoker")
obj_smokers = obj[featuresToKeep,samplesToKeep]
obj_smokers
head(pData(obj_smokers),3)

## ----normFactors------------------------------------------
head(normFactors(obj))
normFactors(obj) <- rnorm(ncol(obj))
head(normFactors(obj))

## ----libSize----------------------------------------------
head(libSize(obj))
libSize(obj) <- rnorm(ncol(obj))
head(libSize(obj))

## ----filterData-------------------------------------------
data(mouseData)
filterData(mouseData,present=10,depth=1000)

## ----mergeMRexperiment------------------------------------
data(mouseData)
newobj = mergeMRexperiments(mouseData,mouseData)
newobj

## ----calculateNormFactors---------------------------------
data(lungData)
p=cumNormStatFast(lungData)

## ----normalizeData----------------------------------------
lungData = cumNorm(lungData,p=p)

## ----wrenchNorm-------------------------------------------
condition = mouseData$diet
mouseData = wrenchNorm(mouseData,condition=condition)

## ----saveData---------------------------------------------
mat = MRcounts(lungData,norm=TRUE,log=TRUE)[1:5,1:5]
exportMat(mat,file=file.path(dataDirectory,"tmp.tsv"))

## ----exportStats------------------------------------------
exportStats(lungData[,1:5],file=file.path(dataDirectory,"tmp.tsv"))
head(read.csv(file=file.path(dataDirectory,"tmp.tsv"),sep="\t"))

## ----removeData, echo=FALSE-------------------------------
system(paste("rm",file.path(dataDirectory,"tmp.tsv")))

## ----fitFeatureModel--------------------------------------
data(lungData)
lungData = lungData[,-which(is.na(pData(lungData)$SmokingStatus))]
lungData=filterData(lungData,present=30,depth=1)
lungData <- cumNorm(lungData, p=.5)
pd <- pData(lungData)
mod <- model.matrix(~1+SmokingStatus, data=pd)
lungres1 = fitFeatureModel(lungData,mod)
head(MRcoefs(lungres1))

## ----preprocess,dev='pdf',out.width='.55\\linewidth',out.height='.55\\linewidth',fig.cap='Relative difference for the median difference in counts from the reference.',fig.align='center',warning=FALSE----
data(lungData)
controls = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-controls]
rareFeatures = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-rareFeatures,]
lungp = cumNormStat(lungTrim,pFlag=TRUE,main="Trimmed lung data")
lungTrim = cumNorm(lungTrim,p=lungp)

## ----zigTesting-------------------------------------------
smokingStatus = pData(lungTrim)$SmokingStatus
bodySite = pData(lungTrim)$SampleType
normFactor = normFactors(lungTrim)
normFactor = log2(normFactor/median(normFactor) + 1)
mod = model.matrix(~smokingStatus+bodySite + normFactor)
settings = zigControl(maxit=10,verbose=TRUE)
fit = fitZig(obj = lungTrim,mod=mod,useCSSoffset = FALSE, 
             control=settings)

# The default, useCSSoffset = TRUE, automatically includes the CSS scaling normalization factor.

## ----contrasts--------------------------------------------
# maxit=1 is for demonstration purposes
settings = zigControl(maxit=1,verbose=FALSE)
mod = model.matrix(~bodySite)
colnames(mod) = levels(bodySite)
# fitting the ZIG model
res = fitZig(obj = lungTrim,mod=mod,control=settings)
# The output of fitZig contains a list of various useful items. hint: names(res).
# 
# Probably the most useful is the limma 'MLArrayLM' object called fit.
zigFit = slot(res,"fit")
finalMod = slot(res,"fit")$design

contrast.matrix = makeContrasts(BAL.A-BAL.B,OW-PSB,levels=finalMod)
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 = eBayes(fit2)
topTable(fit2)

# See help pages on decideTests, topTable, topTableF, vennDiagram, etc.

## ----fittedResult,tidy=TRUE-------------------------------
taxa = 
  sapply(strsplit(as.character(fData(lungTrim)$taxa),split=";"),
         function(i){i[length(i)]})
head(MRcoefs(fit,taxa=taxa,coef=2))

## ----timeSeries-------------------------------------------
# vignette("fitTimeSeries")

## ----perm-------------------------------------------------
coeffOfInterest = 2
res = fitLogNormal(obj = lungTrim, mod = mod, useCSSoffset = FALSE, B = 10, coef = coeffOfInterest)

# extract p.values and adjust for multiple testing
# res$p are the p-values calculated through permutation
adjustedPvalues = p.adjust(res$p,method="fdr")

# extract the absolute fold-change estimates
foldChange = abs(res$fit$coef[,coeffOfInterest])

# determine features still significant and order by the 
sigList = which(adjustedPvalues <= .05)
sigList = sigList[order(foldChange[sigList])]

# view the top taxa associated with the coefficient of interest.
head(taxa[sigList])

## ----presenceAbsence--------------------------------------
classes = pData(mouseData)$diet
res = fitPA(mouseData[1:5,],cl=classes)
# Warning - the p-value is calculating 1 despite a high odd's ratio.
head(res)

## ----discOdds---------------------------------------------
classes = pData(mouseData)$diet
res = fitDO(mouseData[1:100,],cl=classes,norm=FALSE,log=FALSE)
head(res)

## ----corTest----------------------------------------------
cors = correlationTest(mouseData[55:60,],norm=FALSE,log=FALSE)
head(cors)

## ----uniqueFeatures---------------------------------------
cl = pData(mouseData)[["diet"]]
uniqueFeatures(mouseData,cl,nsamples = 10,nreads = 100)

## ----aggTax-----------------------------------------------
obj = aggTax(mouseData,lvl='phylum',out='matrix')
head(obj[1:5,1:5])

## ----aggSamp----------------------------------------------
obj = aggSamp(mouseData,fct='mouseID',out='matrix')
head(obj[1:5,1:5])

## ----interactiveDisplay-----------------------------------
# Calling display on the MRexperiment object will start a browser session with interactive plots.

# require(interactiveDisplay)
# display(mouseData)

## ----heatmapData,fig.cap='Left) Abundance heatmap (plotMRheatmap). Right) Correlation heatmap (plotCorr).',dev='pdf',fig.show='hold',out.width='.5\\linewidth', out.height='.5\\linewidth'----
trials = pData(mouseData)$diet
heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)

# plotMRheatmap
plotMRheatmap(obj=mouseData,n=200,cexRow = 0.4,cexCol = 0.4,trace="none",
                col = heatmapCols,ColSideColors = heatmapColColors)

# plotCorr
plotCorr(obj=mouseData,n=200,cexRow = 0.25,cexCol = 0.25,
         trace="none",dendrogram="none",col=heatmapCols)

## ----MDSandRareplots,fig.cap='Left) CMDS of features (plotOrd). Right) Rarefaction effect (plotRare).',dev='pdf',fig.show='hold',out.width='.5\\linewidth', out.height='.5\\linewidth'----
cl = factor(pData(mouseData)$diet)

# plotOrd - can load vegan and set distfun = vegdist and use dist.method="bray"
plotOrd(mouseData,tran=TRUE,usePCA=FALSE,useDist=TRUE,bg=cl,pch=21)

# plotRare
res = plotRare(mouseData,cl=cl,pch=21,bg=cl)

# Linear fits for plotRare / legend
tmp=lapply(levels(cl), function(lv) 
  lm(res[,"ident"]~res[,"libSize"]-1, subset=cl==lv))
for(i in 1:length(levels(cl))){
   abline(tmp[[i]], col=i)
}
legend("topleft", c("Diet 1","Diet 2"), text.col=c(1,2),box.col=NA)

## ----plotOTUData,fig.cap='Left) Abundance plot (plotOTU). Right) Multiple OTU abundances (plotGenus).',dev='pdf',fig.show='hold',out.width='.5\\linewidth', out.height='.5\\linewidth',tidy=TRUE----
head(MRtable(fit,coef=2,taxa=1:length(fData(lungTrim)$taxa)))
patients=sapply(strsplit(rownames(pData(lungTrim)),split="_"),
          function(i){
            i[3]
          })
pData(lungTrim)$patients=patients
classIndex=list(smoker=which(pData(lungTrim)$SmokingStatus=="Smoker"))
classIndex$nonsmoker=which(pData(lungTrim)$SmokingStatus=="NonSmoker")
otu = 779

# plotOTU
plotOTU(lungTrim,otu=otu,classIndex,main="Neisseria meningitidis")

# Now multiple OTUs annotated similarly
x = fData(lungTrim)$taxa[otu]
otulist = grep(x,fData(lungTrim)$taxa)

# plotGenus
plotGenus(lungTrim,otulist,classIndex,labs=FALSE,
            main="Neisseria meningitidis")

lablist<- c("S","NS")
axis(1, at=seq(1,6,by=1), labels = rep(lablist,times=3))

## ----plotFeatureData,fig.cap='Plot of raw abundances',dev='pdf',fig.show='hold',out.width='.5\\linewidth', out.height='.5\\linewidth',tidy=TRUE----
classIndex=list(Western=which(pData(mouseData)$diet=="Western"))
classIndex$BK=which(pData(mouseData)$diet=="BK")
otuIndex = 8770

# par(mfrow=c(1,2))
dates = pData(mouseData)$date
plotFeature(mouseData,norm=FALSE,log=FALSE,otuIndex,classIndex,
            col=dates,sortby=dates,ylab="Raw reads")

## ----cite-------------------------------------------------
citation("metagenomeSeq")

## ----sessionInfo------------------------------------------
sessionInfo()

