context("testing the effect of seed on toy example")
seed_1 <- 1
seed_2 <- 2
# test two seeds are different
expect_false(isTRUE(equals(seed_1, seed_2)))
set.seed(seed_1)
res_1 <- sample(LETTERS, 5)

set.seed(seed_2)
res_2 <- sample(LETTERS, 5)

set.seed(seed_1)
res_3 <- sample(LETTERS, 5)

# test using two different seeds: we expect two different results
expect_false(isTRUE(all.equal(res_1, res_2)))

# test using same seeds: we expect same results
expect_true(isTRUE(all.equal(res_1, res_3)))

context("testing the effect of seed on somnibus functions")

data(RAdat)
# this should be done within BSMethEM
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
ref_mut$est.pi[1000] <- -999.0

set.seed(seed_1)
out_1 <- BSMethEM(data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon, epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail)

set.seed(seed_2)
out_2 <- BSMethEM(data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon, epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail)

set.seed(seed_1)
out_3 <- BSMethEM(data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon, epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail)

# test using two different seeds: we expect two different results
# but because seeds don't affect BSMethEM, and later mgcv function called we are expecting the same result
expect_true(isTRUE(all.equal(out_1, out_2)))

# test using same seeds: we expect same results
expect_true(isTRUE(all.equal(out_1, out_3)))


context("testing gam function from mgcv package")
# example taken from https://cran.r-project.org/web/packages/mgcv/mgcv.pdf
set.seed(seed_1)
dat <- gamSim(5, n = 10, scale = 2)

b_1 <- gam(y ~ x0 + s(x1) + s(x2), data = dat, seed = seed_1)

b_2 <- gam(y ~ x0 + s(x1) + s(x2), data = dat, seed = seed_2)

b_3 <- gam(y ~ x0 + s(x1) + s(x2), data = dat, seed = seed_1)

# test using two different seeds: we expect two different results
# but because seeds don't affect mgcv functions we are expecting the same result
expect_true(isTRUE(all.equal(b_1$linear.predictors, b_2$linear.predictors)))

# test using same seeds: we expect same results
expect_true(isTRUE(all.equal(b_1$linear.predictors, b_3$linear.predictors)))
