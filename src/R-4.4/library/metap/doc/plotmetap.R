### R code from vignette source 'plotmetap.Rnw'

###################################################
### code chunk number 1: plotmetap.Rnw:70-76
###################################################
library(metap)
data(dat.metap)
teach <- dat.metap$teachexpect
validity <- dat.metap$validity$p
zhang <- dat.metap$zhang
print(validity)


###################################################
### code chunk number 2: plotp
###################################################
plotp(validity, main = "Validity data")


###################################################
### code chunk number 3: plotfunc
###################################################
plotp(validity, main = "Validity data", plotversion = "old")


###################################################
### code chunk number 4: teachlinear
###################################################
plotp(teach)


###################################################
### code chunk number 5: teachlog
###################################################
plotp(teach, log10 = TRUE)


###################################################
### code chunk number 6: plotmetap.Rnw:167-168
###################################################
logitp(teach)


###################################################
### code chunk number 7: simple
###################################################
schweder(validity)


###################################################
### code chunk number 8: withlines
###################################################
schweder(validity, drawline = c("bh", "ls", "ab"),
   ls.control = list(frac = 0.5), ab.control = list(a = 0, b = 0.01))


###################################################
### code chunk number 9: albatros
###################################################
validity <- dat.metap$validity
fit.v <- albatros(validity$p, validity$n,
   contours = list(type = "corr", contvals = c(0.25, 0.5, 0.8), ltys = 1:3),
      axes = list(ylimit = c(1,200),  lefttext = "Negative correlation",
         righttext = "Positive correlation"),
   main = "Validity")


###################################################
### code chunk number 10: zhang
###################################################
fit.z <- albatros(zhang$p, zhang$n,
   contours = list(type = "smd", contvals = c(0.25, 0.5, 1), ltys = 1:3),
   plotpars = list(pchs = letters[unclass(dat.metap$zhang$phase)]),
   axes = list(lefttext = "Favours control", righttext = "Favours exercise"),
   main = "Zhang"
   )


