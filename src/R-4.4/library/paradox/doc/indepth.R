## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment = "#>", collapse = TRUE)

## -----------------------------------------------------------------------------
library("paradox")

ps1 = ps(a = p_int(init = 1))
ps2 = ps1
ps3 = ps1$clone(deep = TRUE)
print(ps1) # the same for ps2 and ps3

## -----------------------------------------------------------------------------
ps1$values$a = 2

## -----------------------------------------------------------------------------
print(ps1) # ps1 value of 'a' was changed
print(ps2) # contains the same reference as ps1, so also changed
print(ps3) # is a "clone" of the old ps1 with 'a' == 1

## -----------------------------------------------------------------------------
library("paradox")
param_set = ps(
  parA = p_lgl(init = FALSE),
  parB = p_int(lower = 0, upper = 10, tags = c("tag1", "tag2")),
  parC = p_dbl(lower = 0, upper = 4, special_vals = list(NULL)),
  parD = p_fct(levels = c("x", "y", "z"), default = "y"),
  parE = p_uty(custom_check = function(x) checkmate::checkFunction(x))
)
param_set

## -----------------------------------------------------------------------------
param_set$lower
param_set$parD$levels
param_set$class

## -----------------------------------------------------------------------------
as.data.table(param_set)

## -----------------------------------------------------------------------------
param_set$test(list(parA = FALSE, parB = 0))
param_set$test(list(parA = "FALSE"))
param_set$check(list(parA = "FALSE"))

## -----------------------------------------------------------------------------
ps1 = ParamSet$new(list(x = p_int(), y = p_dbl()))
ps2 = ParamSet$new(list(z = p_fct(levels = c("a", "b", "c"))))
ps_all = c(ps1, ps2)
print(ps_all)
ps_all$subset(c("x", "z"))

## -----------------------------------------------------------------------------
as.data.table(ps_all)

## -----------------------------------------------------------------------------
ps1$values = list(x = 1, y = 1.5)
ps1$values$y = 2.5
print(ps1$values)

## ----error = TRUE-------------------------------------------------------------
ps1$values$x = 1.5

## -----------------------------------------------------------------------------
p = ps(
  A = p_lgl(init = FALSE),
  B = p_int(lower = 0, upper = 10, depends = D %in% c("x", "y")),
  C = p_dbl(lower = 0, upper = 4),
  D = p_fct(levels = c("x", "y", "z"), depends = A == FALSE)
)

## -----------------------------------------------------------------------------
p$check(list(A = FALSE, D = "x", B = 1), check_strict = TRUE)  # OK: all dependencies met
p$check(list(A = FALSE, D = "z", B = 1), check_strict = TRUE)  # B's dependency is not met
p$check(list(A = FALSE, B = 1), check_strict = TRUE)           # B's dependency is not met
p$check(list(A = FALSE, D = "z"), check_strict = TRUE)         # OK: B is absent
p$check(list(A = TRUE), check_strict = TRUE)                   # OK: neither B nor D present
p$check(list(A = TRUE, D = "x", B = 1), check_strict = TRUE)   # D's dependency is not met
p$check(list(A = TRUE, B = 1), check_strict = TRUE)            # B's dependency is not met

## -----------------------------------------------------------------------------
p$deps

## -----------------------------------------------------------------------------
ps2d = ps_replicate(ps(x = p_dbl(lower = 0, upper = 1)), 2)
print(ps2d)

## -----------------------------------------------------------------------------
ps2d = ps_replicate(ps(x = p_dbl(0, 1), y = p_int(0, 10)), 2, tag_params = TRUE)
ps2d$values = list(rep1.x = 0.2, rep2.x = 0.4, rep1.y = 3, rep2.y = 4)
ps2d$tags
ps2d$get_values(tags = "param_x")

## -----------------------------------------------------------------------------
ps_small = ps(A = p_dbl(0, 1), B = p_dbl(0, 1))
design = generate_design_grid(ps_small, 2)
print(design)

## -----------------------------------------------------------------------------
generate_design_grid(ps_small, param_resolutions = c(A = 3, B = 2))

## -----------------------------------------------------------------------------
pvrand = generate_design_random(ps_small, 500)
pvlhs = generate_design_lhs(ps_small, 500)
pvsobol = generate_design_sobol(ps_small, 500)

## ----echo = FALSE, out.width="45%", fig.show = "hold", fig.width = 4, fig.height = 4----
par(mar=c(4, 4, 2, 1))
plot(pvrand$data, main = "'random' design", xlim = c(0, 1), ylim=c(0, 1))
plot(pvlhs$data, main = "'lhs' design", xlim = c(0, 1), ylim=c(0, 1))
plot(pvsobol$data, main = "'sobol' design", xlim = c(0, 1), ylim=c(0, 1))

## -----------------------------------------------------------------------------
sampA = Sampler1DCateg$new(ps(x = p_fct(letters)))
sampA$sample(5)

## -----------------------------------------------------------------------------
p = ps(
  A = p_lgl(),
  B = p_int(0, 10, depends = A == TRUE)
)

p_subspaces = p$subspaces()

sampH = SamplerHierarchical$new(p,
  list(Sampler1DCateg$new(p_subspaces$A),
    Sampler1DUnif$new(p_subspaces$B))
)
sampled = sampH$sample(1000)
head(sampled$data)
table(sampled$data[, c("A", "B")], useNA = "ifany")

## -----------------------------------------------------------------------------
sampJ = SamplerJointIndep$new(
  list(Sampler1DUnif$new(ps(x = p_dbl(0, 1))),
    Sampler1DUnif$new(ps(y = p_dbl(0, 1))))
)
sampJ$sample(5)

## -----------------------------------------------------------------------------
psexp = ps(par = p_dbl(0, 1, trafo = function(x) -log(x)))

design = generate_design_random(psexp, 3)
print(design)  # not transformed: between 0 and 1
design$transpose()  # trafo is TRUE

## -----------------------------------------------------------------------------
design$transpose(trafo = FALSE)

## -----------------------------------------------------------------------------
psexp = ps(par = p_dbl(0, 1))
psexp$extra_trafo = function(x, param_set) {
  x$par = -log(x$par)
  x
}

## -----------------------------------------------------------------------------
methodPS = ps(fun = p_uty(custom_check = function(x) checkmate::checkFunction(x, nargs = 1)))

print(methodPS)

## -----------------------------------------------------------------------------
samplingPS = ps(
  fun = p_fct(c("mean", "median", "min", "max"),
    trafo = function(x) get(x, mode = "function"))
)

## -----------------------------------------------------------------------------
design = generate_design_random(samplingPS, 2)
print(design)

## -----------------------------------------------------------------------------
xvals = design$transpose()
print(xvals[[1]])

## -----------------------------------------------------------------------------
methodPS$check(xvals[[1]])
xvals[[1]]$fun(1:10)

## -----------------------------------------------------------------------------
samplingPS = ps(
  fun = p_fct(list("mean" = mean, "median" = median, "min" = min, "max" = max))
)

generate_design_random(samplingPS, 1)$transpose()

## -----------------------------------------------------------------------------
samplingPS2 = ps(quantile = p_dbl(0, 1),
  .extra_trafo = function(x, param_set) {
    # x$quantile is a `numeric(1)` between 0 and 1.
    # We want to turn it into a function!
    list(fun = function(input) quantile(input, x$quantile))
  }
)

## -----------------------------------------------------------------------------
design = generate_design_random(samplingPS2, 2)
print(design)

## -----------------------------------------------------------------------------
xvals = design$transpose()
print(xvals[[1]])
methodPS$check(xvals[[1]])
xvals[[1]]$fun(1:10)

## -----------------------------------------------------------------------------
search_space = ps()
print(search_space)

## -----------------------------------------------------------------------------
search_space = ps(
  cost = p_dbl(lower = 0.1, upper = 10),
  kernel = p_fct(levels = c("polynomial", "radial"))
)
print(search_space)

## -----------------------------------------------------------------------------
search_space = ps(cost = p_dbl(0.1, 10), kernel = p_fct(c("polynomial", "radial")))

## -----------------------------------------------------------------------------
library("data.table")
rbindlist(generate_design_grid(search_space, 3)$transpose())

## -----------------------------------------------------------------------------
search_space = ps(
  cost = p_dbl(-1, 1, trafo = function(x) 10^x),
  kernel = p_fct(c("polynomial", "radial"))
)
rbindlist(generate_design_grid(search_space, 3)$transpose())

## -----------------------------------------------------------------------------
search_space = ps(
  cost = p_dbl(-1, 1, trafo = function(x) 10^x),
  kernel = p_fct(c("polynomial", "radial")),
  .extra_trafo = function(x, param_set) {
    if (x$kernel == "polynomial") {
      x$cost = x$cost * 2
    }
    x
  }
)
rbindlist(generate_design_grid(search_space, 3)$transpose())

## -----------------------------------------------------------------------------
search_space = ps(
  class.weights = p_dbl(0.1, 0.9, trafo = function(x) c(spam = x, nonspam = 1 - x))
)
generate_design_grid(search_space, 3)$transpose()

## -----------------------------------------------------------------------------
search_space = ps(
  cost = p_fct(c(0.1, 3, 10)),
  kernel = p_fct(c("polynomial", "radial"))
)
rbindlist(generate_design_grid(search_space, 3)$transpose())

## -----------------------------------------------------------------------------
search_space = ps(
  cost = p_fct(c("0.1", "3", "10"),
    trafo = function(x) list(`0.1` = 0.1, `3` = 3, `10` = 10)[[x]]),
  kernel = p_fct(c("polynomial", "radial"))
)
rbindlist(generate_design_grid(search_space, 3)$transpose())

## -----------------------------------------------------------------------------
search_space = ps(
  cost = p_fct(c(0.1, 3, 10)),
  kernel = p_fct(c("polynomial", "radial"))
)
typeof(search_space$params$cost$levels)

## -----------------------------------------------------------------------------
search_space = ps(
  class.weights = p_fct(
    list(
      candidate_a = c(spam = 0.5, nonspam = 0.5),
      candidate_b = c(spam = 0.3, nonspam = 0.7)
    )
  )
)
generate_design_grid(search_space)$transpose()

## -----------------------------------------------------------------------------
search_space = ps(
  cost = p_dbl(-1, 1, trafo = function(x) 10^x),
  kernel = p_fct(c("polynomial", "radial")),
  degree = p_int(1, 3, depends = kernel == "polynomial")
)
rbindlist(generate_design_grid(search_space, 3)$transpose(), fill = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  library("mlr3learners")
#  learner = lrn("classif.svm")
#  learner$param_set$values$kernel = "polynomial" # for example
#  learner$param_set$values$degree = to_tune(lower = 1, upper = 3)
#  
#  print(learner$param_set$search_space())
#  
#  rbindlist(generate_design_grid(
#    learner$param_set$search_space(), 3)$transpose()
#  )

## ----eval = FALSE-------------------------------------------------------------
#  learner$param_set$values$shrinking = to_tune()
#  
#  print(learner$param_set$search_space())
#  
#  rbindlist(generate_design_grid(
#    learner$param_set$search_space(), 3)$transpose()
#  )

## ----eval = FALSE-------------------------------------------------------------
#  learner$param_set$values$type = "C-classification" # needs to be set because of a bug in paradox
#  learner$param_set$values$cost = to_tune(c(val1 = 0.3, val2 = 0.7))
#  learner$param_set$values$shrinking = to_tune(p_lgl(depends = cost == "val2"))
#  
#  print(learner$param_set$search_space())
#  
#  rbindlist(generate_design_grid(learner$param_set$search_space(), 3)$transpose(), fill = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  learner$param_set$values$cost = NULL
#  learner$param_set$values$shrinking = NULL
#  learner$param_set$values$kernel = to_tune(c("polynomial", "radial"))
#  
#  print(learner$param_set$search_space())
#  
#  rbindlist(generate_design_grid(learner$param_set$search_space(), 3)$transpose(), fill = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  learner$param_set$values$class.weights = to_tune(
#    ps(spam = p_dbl(0.1, 0.9), nonspam = p_dbl(0.1, 0.9),
#      .extra_trafo = function(x, param_set) list(c(spam = x$spam, nonspam = x$nonspam))
#  ))
#  head(generate_design_grid(learner$param_set$search_space(), 3)$transpose(), 3)

