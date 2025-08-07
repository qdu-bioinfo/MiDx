### R code from vignette source 'compare.Rnw'

###################################################
### code chunk number 1: compare.Rnw:178-179
###################################################
cancel <- c(0.001, 0.001, 0.999, 0.999)


###################################################
### code chunk number 2: compare.Rnw:183-202
###################################################
library(metap)
data(dat.metap)
validity <- dat.metap$validity$p
genvec <- function(pvals, kvals, fun, name) {
   ps <- length(pvals)
   ks <- length(kvals)
   temp <- matrix(-1, nrow = ps, ncol = ks)
   for(i in 1:ps)
   for(j in 1:ks) {
      temp[i, j] <- fun(rep(pvals[i], kvals[j]))$p
   }
   temp2 <- as.vector(temp)
   res <- data.frame(method = rep(name, length(temp2)),
      p = rep(pvals, ks),
      k = rep(kvals, each = ps),
      g = temp2
   )
   res
}


###################################################
### code chunk number 3: compare.Rnw:229-237
###################################################
   kvals <- c(4, 5, 6, 8, 10, 15, 20)
   pvals <- c(0.2, 0.3, 0.3679, 0.4, 0.5, 0.6)
   dat <- rbind(
      genvec(pvals, kvals, logitp, "logitp"),
      genvec(pvals, kvals, meanz, "meanz"),
      genvec(pvals, kvals, sumlog, "sumlog"),
      genvec(pvals, kvals, sumz, "sumz")
   )


###################################################
### code chunk number 4: transeqp
###################################################
   lattice::xyplot(g ~ k | method, groups = p, type = "l", data = dat,
      auto.key = list(space = "left", lines = TRUE, title = "p"),
      ylab = "g(p)"
   )


###################################################
### code chunk number 5: compare.Rnw:259-282
###################################################
set.seed(18122019)
temp <- matrix(runif(10000), nrow = 100)
fisher <- apply(temp, 1, function(x) sumlog(x)$p)
lanc4 <- apply(temp, 1, function(x) invchisq(x, 4)$p)
lanc16 <- apply(temp, 1, function(x) invchisq(x, 16)$p)
lanc256 <- apply(temp, 1, function(x) invchisq(x, 256)$p)
banda <- function(x, y) {
   res <- data.frame(sum = x + y, diff = (x - y))
   res
}
dat <- data.frame(rbind(banda(fisher, lanc4),
   banda(fisher, lanc16),
   banda(fisher, lanc256),
   banda(lanc4, lanc16),
   banda(lanc4, lanc256),
   banda(lanc16, lanc256)
   ),
   name = factor(c(rep("FL4", 100), rep("FL16", 100),
      rep("FL256", 100), rep("L4L16", 100),
      rep("L4L256", 100), rep("L16L256", 100)),
      levels = c("FL4", "FL16", "FL256", "L4L16", "L4L256", "L16L256")
   )
)


###################################################
### code chunk number 6: fishlanc
###################################################
   lattice::xyplot(diff ~ sum | name, data = dat,
   panel = function(x, y, ...) {
      lattice::panel.xyplot(x, y, ...)
      lattice::panel.abline(h = mean(y), lty = 2)
      lattice::panel.abline(h = mean(y) + 1.96 * sd(y), lty = 3)
      lattice::panel.abline(h = mean(y) - 1.96 * sd(y), lty = 3)
   }
   )


###################################################
### code chunk number 7: compare.Rnw:296-317
###################################################
stouff <- apply(temp, 1, function(x) sumz(x)$p)
invt4 <- apply(temp, 1, function(x) invt(x, 4)$p)
invt16 <- apply(temp, 1, function(x) invt(x, 16)$p)
invt256 <- apply(temp, 1, function(x) invt(x, 256)$p)
banda <- function(x, y) {
   res <- data.frame(sum = x + y, diff = (x - y))
   res
}
dat <- data.frame(rbind(banda(stouff, invt4),
   banda(stouff, invt16),
   banda(stouff, invt256),
   banda(invt4, invt16),
   banda(invt4, invt256),
   banda(invt16, invt256)
   ),
   name = factor(c(rep("St4", 100), rep("St16", 100),
      rep("St256", 100), rep("t4t16", 100),
      rep("t4t256", 100), rep("t16t256", 100)),
      levels = c("St4", "St16", "St256", "t4t16", "t4t256", "t16t256")
   )
)


###################################################
### code chunk number 8: stouffinvt
###################################################
lattice::xyplot(diff ~ sum | name, data = dat,
   panel = function(x, y, ...) {
      lattice::panel.xyplot(x, y, ...)
      lattice::panel.abline(h = mean(y), lty = 2)
      lattice::panel.abline(h = mean(y) + 1.96 * sd(y), lty = 3)
      lattice::panel.abline(h = mean(y) - 1.96 * sd(y), lty = 3)
   }
)


###################################################
### code chunk number 9: compare.Rnw:400-402
###################################################
meanz(c(0.3, 0.31))$p
meanz(c(0.1, 0.2))$p


###################################################
### code chunk number 10: compare.Rnw:438-442
###################################################
log10p <- function(x) {
   res <- round(-log(x, base = 10), 2)
   res
}


###################################################
### code chunk number 11: compare.Rnw:492-501
###################################################
   kvals <- c(4, 5, 6, 8, 10, 15, 20)
   pvals <- c(0.2, 0.3, 0.3679, 0.4, 0.5, 0.6)
   dat <- rbind(
      genvec(pvals, kvals, meanp, "meanp"),
      genvec(pvals, kvals, maximump, "maximump"),
      genvec(pvals, kvals, minimump, "minimump"),
      genvec(pvals, kvals, sump, "sump"),
      genvec(pvals, kvals, votep, "votep")
   )


###################################################
### code chunk number 12: untranseqp
###################################################
   lattice::xyplot(g ~ k | method, groups = p, type = "l", data = dat,
      auto.key = list(space = "left", lines = TRUE, title = "p"),
      ylab = "g(p)"
   )


