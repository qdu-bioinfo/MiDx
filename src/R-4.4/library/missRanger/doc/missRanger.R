## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(missRanger)

set.seed(3)

iris_NA <- generateNA(iris, p = 0.1)
head(iris_NA)
 
imp <- missRanger(iris_NA, num.trees = 100)
head(imp)

## -----------------------------------------------------------------------------
imp <- missRanger(iris_NA, pmm.k = 5, num.trees = 100, verbose = 0)
head(imp)

## -----------------------------------------------------------------------------
imp <- missRanger(iris_NA, pmm.k = 5, num.trees = 200, mtry = 1, verbose = 0)

## -----------------------------------------------------------------------------
imp <- missRanger(
  iris_NA, pmm.k = 5, num.trees = 100, keep_forests = TRUE, verbose = 0
)
imp

summary(imp)

# Out-of-sample application
# saveRDS(imp, file = "imputation_model.rds")
# imp <- readRDS("imputation_model.rds")
predict(imp, head(iris_NA))

## -----------------------------------------------------------------------------
# Impute all variables with all (default)
m <- missRanger(iris_NA, formula = . ~ ., pmm.k = 5, num.trees = 100, verbose = 0)

# Don't use Species for imputation
m <- missRanger(iris_NA, . ~ . - Species, pmm.k = 5, num.trees = 100, verbose = 0)

# Impute Sepal.Length by Species (or not?)
m <- missRanger(iris_NA, Sepal.Length ~ Species, pmm.k = 5, num.trees = 100)
head(m)

# Only univariate imputation was done! Why? Because Species contains missing values
# itself and needs to appear on the LHS as well:
m <- missRanger(iris_NA, Sepal.Length + Species ~ Species, pmm.k = 5, num.trees = 100)
head(m)

# Impute all variables univariately
m <- missRanger(iris_NA, . ~ 1, verbose = 0)

