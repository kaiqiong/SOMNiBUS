context("testing the effect of seed on toy example")

# the only place that requires seed is in the function extractDesignMatrix

#  sub.samp <- max(1000, 2 * length(GamObj$coefficients))
#if (nrow(GamObj$model) > sub.samp) {
    ## subsample to get X for p-values calc.  sample these rows from X
#   ind <- sample(seq_len(nrow(GamObj$model)), sub.samp, replace=FALSE)
#    X_d <- predict(GamObj, GamObj$model[ind, ], type="lpmatrix")
#}

# the sample line only got executed when the number of free parameters are greater than 1000,
# which is not the case in the example of RAdat, where only 58 (sum(n.k) + n_samples)
# parameters are present. So seed shouldn't influence the results for RAdat

library(mgcv)

seed_1 <- 1
seed_2 <- 2

test_that("the two seeds are different", {
  expect_false(isTRUE(equals(seed_1, seed_2)))
})

set.seed(seed_1)
res_1 <- sample(LETTERS, 5)

set.seed(seed_2)
res_2 <- sample(LETTERS, 5)

set.seed(seed_1)
res_3 <- sample(LETTERS, 5)

test_that("using two different seeds we get two different results", {
  expect_false(isTRUE(all.equal(res_1, res_2)))
})

test_that("using the same seeds we get the same results", {
  expect_true(isTRUE(all.equal(res_1, res_3)))
})

context("testing the effect of seed on somnibus functions")

data(RAdat)
# this should be done within binomRegMethModel
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])

# default parameters to be used
epsilon <- 10^(-6)
epsilon.lambda <- 10^(-3)
p0 <- 0.003307034
p1 <- 0.9
maxStep <- 200
detail <- FALSE
n.k <- rep(5, 3)

path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_results.RDS", sep = "")
ref <- readRDS(path_ref_data)
ref_mut <- readRDS(path_ref_data)
ref_mut$est.pi[1000] <- -999

set.seed(seed_1)
out_1 <- binomRegMethModel(
  data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon,
  epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail
)

set.seed(seed_2)
out_2 <- binomRegMethModel(
  data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon,
  epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail
)

set.seed(seed_1)
out_3 <- binomRegMethModel(
  data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon,
  epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail
)

test_that("using two different seeds we get the same results with binomRegMethModel", {
  # this is wrong different seeds settings should output different results with
  # binomRegMethModel
  expect_true(isTRUE(all.equal(out_1, out_2)))
})

test_that("using the same seed we get the same results with binomRegMethModel", {
  # test using same seeds: we expect same results
  expect_true(isTRUE(all.equal(out_1, out_3)))
})

context("testing gam function from mgcv package")
# example taken from https://cran.r-project.org/web/packages/mgcv/mgcv.pdf
set.seed(seed_1)
dat <- gamSim(5, n = 10, scale = 2)

b_1 <- gam(y ~ x0 + s(x1) + s(x2), data = dat, seed = seed_1)

b_2 <- gam(y ~ x0 + s(x1) + s(x2), data = dat, seed = seed_2)

b_3 <- gam(y ~ x0 + s(x1) + s(x2), data = dat, seed = seed_1)

test_that("using two different seeds we get the same results with mgcv::gam", {
  expect_true(isTRUE(all.equal(b_1$linear.predictors, b_2$linear.predictors)))
})

test_that("using same seeds we get same results using with mgcv::gam", {
  expect_true(isTRUE(all.equal(b_1$linear.predictors, b_3$linear.predictors)))
})
