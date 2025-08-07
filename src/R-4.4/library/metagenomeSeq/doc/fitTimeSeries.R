## ----include=FALSE------------------------------------------------------------
require(knitr)
opts_chunk$set(concordance=TRUE,tidy=TRUE)

## ----config,echo=FALSE-----------------------------------------
options(width = 65)
options(continue=" ")
options(warn=-1)
set.seed(42)

## ----requireMetagenomeSeq,warning=FALSE,message=FALSE----------
library(metagenomeSeq)
library(gss)

## ----dataset2,tidy=FALSE---------------------------------------
data(mouseData)
mouseData

## ----createMRexperiment1---------------------------------------
# Creating mock sample replicates
sampleID = rep(paste("sample",1:10,sep=":"),times=20)
# Creating mock class membership
class = rep(c(rep(0,5),rep(1,5)),times=20)
# Creating mock time
time = rep(1:20,each=10)

phenotypeData = AnnotatedDataFrame(data.frame(sampleID,class,time))
# Creating mock abundances
set.seed(1)
# No difference
measurement1 = rnorm(200,mean=100,sd=1)
# Some difference
measurement2 = rnorm(200,mean=100,sd=1)
measurement2[1:5]=measurement2[1:5] + 100
measurement2[11:15]=measurement2[11:15] + 100
measurement2[21:25]=measurement2[21:25] + 50
mat = rbind(measurement1,measurement2)
colnames(mat) = 1:200
mat[1:2,1:10]

## ----createMRexperiment2---------------------------------------
# This is an example of potential lvl's to aggregate by.
data(mouseData)
colnames(fData(mouseData))

## ----createMRexperiment3,tidy=FALSE----------------------------
obj = newMRexperiment(counts=mat,phenoData=phenotypeData)
obj
res1 = fitTimeSeries(obj,feature=1,
              class='class',time='time',id='sampleID',
              B=10,norm=FALSE,log=FALSE)
res2 = fitTimeSeries(obj,feature=2,
              class='class',time='time',id='sampleID',
              B=10,norm=FALSE,log=FALSE)


classInfo = factor(res1$data$class)

## ----plotMRexperiment3,tidy=FALSE------------------------------
par(mfrow=c(3,1))
plotClassTimeSeries(res1,pch=21,bg=classInfo)
plotTimeSeries(res2)
plotClassTimeSeries(res2,pch=21,bg=classInfo)

## ----timeSeries------------------------------------------------
res = fitTimeSeries(obj=mouseData,lvl="class",feature="Actinobacteria",class="status",id="mouseID",time="relativeTime",B=10)

# We observe a time period of differential abundance for "Actinobacteria"
res$timeIntervals

str(res)

## ----timeSeriesAllClasses, tidy=FALSE--------------------------
set.seed(123)
classes = unique(fData(mouseData)[,"class"])

timeSeriesFits = lapply(classes,function(i){
        fitTimeSeries(obj=mouseData,
            feature=i,
            class="status",
            id="mouseID",
            time="relativeTime",
            lvl='class',
            C=.3,# a cutoff for 'interesting' 
            B=1) # B is the number of permutations and should clearly not be 1
    })
names(timeSeriesFits) = classes

# Removing classes of bacteria without a potentially
# interesting time interval difference.
timeSeriesFits = lapply(timeSeriesFits,function(i){i[[1]]})[-grep("No",timeSeriesFits)]

# Naming the various interesting time intervals.
for(i in 1:length(timeSeriesFits)){
    rownames(timeSeriesFits[[i]]) = 
      paste(
        paste(names(timeSeriesFits)[i]," interval",sep=""),
        1:nrow(timeSeriesFits[[i]]),sep=":"
      )
}

# Merging into a table.
timeSeriesFits = do.call(rbind,timeSeriesFits)

# Correcting for multiple testing.
pvalues = timeSeriesFits[,"p.value"]
adjPvalues = p.adjust(pvalues,"bonferroni")
timeSeriesFits = cbind(timeSeriesFits,adjPvalues)

head(timeSeriesFits)

## ----timeSeriesPlotting----------------------------------------
par(mfrow=c(2,1))
plotClassTimeSeries(res,pch=21,
                    bg=res$data$class,ylim=c(0,8))
plotTimeSeries(res)

## ----cite------------------------------------------------------
citation("metagenomeSeq")

## ----sessionInfo-----------------------------------------------
sessionInfo()

