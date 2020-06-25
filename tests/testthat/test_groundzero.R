require(mgcv)
require(SOMNiBUS)

context("testing that we recover the same values over and over with a given dataset")

seed_1 <- 1

data(RAdat)
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

set.seed(seed_1)
out_1 <- BSMethEM(data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon, epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail)
# saveRDS(out_1, path_ref_data)

# test all values in ref equals actual values coming from BSMethEM because both result were obtained using the same seed
expect_true(isTRUE(all.equal(out_1, ref)))
