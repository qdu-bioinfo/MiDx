### R code from vignette source 'metap.Rnw'

###################################################
### code chunk number 1: metap.Rnw:119-120
###################################################
library(metap)


###################################################
### code chunk number 2: metap.Rnw:137-141
###################################################
pvals <- c(0.1, 0.1, 0.9, 0.9, 0.9, 0.9)
istwo <- c(TRUE,  FALSE, TRUE, FALSE, TRUE, FALSE)
toinvert <- c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
two2one(pvals, two = istwo, invert = toinvert)


###################################################
### code chunk number 3: plotp
###################################################
data(dat.metap)
validity <- dat.metap$validity$p
plotp(validity)


###################################################
### code chunk number 4: metap.Rnw:182-183
###################################################
print(validity)


###################################################
### code chunk number 5: metap.Rnw:212-213
###################################################
sumlog(validity)


