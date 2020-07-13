context("testing BSMethSim")
fp <- 0.003307034 # False positive rate (float)
tp <- 0.9 # True positive rate (float)

# loading reference results
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_results.RDS", sep = "")

# loading
ref <- readRDS(path_ref_data)
data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# theta, beta parameters vectors are draw from normal distribution (mu = 0, std=1)
# coverage
n_samples <- length(unique(RAdat.f$ID)) # int: number of samples
n_positions <- length(ref$uni.pos) # int: number of CpG positions
genomics_positions <- ref$uni.pos # numeric vector of (int) of n_positions
theta <- matrix(rnorm(n_positions, 0, 1), n_positions, 1) # numeric vector (floats) of size n_positions : GAM parameters vector
beta <- matrix(rnorm(n_positions * 2, 0, 1), n_positions, 2) # numeric vetor (floats) of size n_positions : GAM parameters vector
coverage <- matrix(sample(1:500, n_samples * n_positions, replace = TRUE), n_samples, n_positions) # numeric matrix (floats or integers?) of size n_samples per n_positions: reads coverage
covariates <- matrix(sample(0:1, n_samples * 2, replace = TRUE), n_samples, 2) # numeric matrix (floats) of size n_samples per n_covariates: covariates matrix
error <- TRUE # adding normaly distributed noise to the simulations (bool)
error_mu <- 0.0 # mean value (float)
error_var <- 1.0 # variance value (float)
link_fct <- "logit" # link fonction (str)

out_BSMethSim <- BSMethSim(n = n_samples, posit = genomics_positions, theta.0 = theta, beta = beta, random.eff = error, mu.e = error_mu, sigma.ee = error_var, p0 = fp, p1 = tp, X = coverage, Z = covariates, binom.link = link_fct)

test_that("the output corresponds to BSMethSim method description", {
  expect_false(length(out_BSMethSim) == 3)
  expect_equal(ncol(out_BSMethSim$S), n_positions)
  expect_equal(nrow(out_BSMethSim$S), n_samples)
  expect_equal(ncol(out_BSMethSim$Y), n_positions)
  expect_equal(nrow(out_BSMethSim$Y), n_samples)
  expect_equal(ncol(out_BSMethSim$theta), n_positions)
  expect_equal(nrow(out_BSMethSim$theta), n_samples)
})
