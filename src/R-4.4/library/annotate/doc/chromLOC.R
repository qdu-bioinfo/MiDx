## ----buildCL, message=FALSE---------------------------------------------------
library("annotate")
z <- buildChromLocation("hgu95av2")
z

## ----showBasicMethods---------------------------------------------------------
organism(z)

dataSource(z)

## The chromLocs list is extremely large. Let's only
## look at one of the elements.
names(chromLocs(z))
chromLocs(z)[["Y"]]

get("32972_at", probesToChrom(z))

chromInfo(z)

get("32972_at", geneSymbols(z))

## ----nChrom-------------------------------------------------------------------
nChrom(z)

