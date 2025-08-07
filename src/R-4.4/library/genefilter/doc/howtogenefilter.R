## ----closeg, message = FALSE--------------------------------------------------
library("Biobase") 
library("genefilter")
data(sample.ExpressionSet)
varLabels(sample.ExpressionSet)
table(sample.ExpressionSet$sex)
table(sample.ExpressionSet$type)

## ----message=FALSE------------------------------------------------------------
f1 <- kOverA(5, 200)
ffun <- filterfun(f1)
wh1 <- genefilter(exprs(sample.ExpressionSet), ffun)
sum(wh1)

## -----------------------------------------------------------------------------
f2 <- ttest(sample.ExpressionSet$type, p=0.1)
wh2 <- genefilter(exprs(sample.ExpressionSet), filterfun(f2))
sum(wh2)

## ----gene-indexing------------------------------------------------------------
ffun_combined <- filterfun(f1, f2)
wh3 <- genefilter(exprs(sample.ExpressionSet), ffun_combined)
sum(wh3)

## ----knnCV--------------------------------------------------------------------
knnCV <- function(x, selectfun, cov, Agg, pselect = 0.01, scale=FALSE) {
   nc <- ncol(x)
   outvals <- rep(NA, nc)
   for(i in seq_len(nc)) {
      v1 <- x[,i]
      expr <- x[,-i]
      glist <- selectfun(expr, cov[-i], p=pselect)
      expr <- expr[glist,]
      if( scale ) {
        expr <- scale(expr)
        v1 <- as.vector(scale(v1[glist]))
      }
      else
         v1 <- v1[glist]
      out <- paste("iter ",i, " num genes= ", sum(glist), sep="")
      print(out)
      Aggregate(row.names(expr), Agg)
      if( length(v1) == 1)
         outvals[i] <- knn(expr, v1, cov[-i], k=5)
      else
          outvals[i] <- knn(t(expr), v1, cov[-i], k=5)
    }
    return(outvals)
}

## ----aggregate1---------------------------------------------------------------
gfun <- function(expr, cov, p=0.05) {
   f2 <- ttest(cov, p=p)
   ffun <- filterfun(f2)
   which <- genefilter(expr, ffun)
  }

## ----aggregate2, results="hide"-----------------------------------------------
library("class")
##scale the genes
##genescale is a slightly more flexible "scale"
##work on a subset -- for speed only
geneData <- genescale(exprs(sample.ExpressionSet)[1:75,], 1) 
Agg <- new("aggregator") 
testcase <- knnCV(geneData, gfun, sample.ExpressionSet$type, 
       Agg,pselect=0.05)  

## ----aggregate3---------------------------------------------------------------
sort(sapply(aggenv(Agg), c), decreasing=TRUE)

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

