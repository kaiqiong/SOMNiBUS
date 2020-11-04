context("testing binomRegMethModelSim")
fp <- 0.003307034 # False positive rate (float)
tp <- 0.9 # True positive rate (float)

# loading reference results
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_results.RDS", sep = "")

# loading
ref <- readRDS(path_ref_data)
data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# theta, beta parameters vectors are draw from normal distribution (mu=0,
# std=1) coverage
n_samples <- length(unique(RAdat.f$ID)) # int: number of samples
n_positions <- length(ref$uni.pos) # int: number of CpG positions
genomics_positions <- ref$uni.pos # numeric vector of (int) of n_positions
theta <- matrix(rnorm(n_positions, 0, 1), n_positions, 1) # numeric vector (floats) of size n_positions : GAM parameters vector
beta <- matrix(rnorm(n_positions * 2, 0, 1), n_positions, 2) # numeric vetor (floats) of size n_positions : GAM parameters vector
coverage <- matrix(
  sample(1:500, n_samples * n_positions, replace = TRUE), n_samples,
  n_positions
)+4 # numeric matrix (floats or integers?) of size n_samples per n_positions: reads coverage
covariates <- matrix(sample(0:1, n_samples * 2, replace = TRUE), n_samples, 2) # numeric matrix (floats) of size n_samples per n_covariates: covariates matrix
error <- TRUE # adding normaly distributed noise (or subject-level random Effect) to the simulations (bool)
error_mu <- 0 # mean value (float)
error_var <- 1 # variance value (float)
link_fct <- "logit" # link fonction (str)
phi_v <- rep(ref$phi_fletcher, n_positions) # numeric vector (floats) of size n_positions:

out_binomRegMethModelSim <- binomRegMethModelSim(
  n = n_samples, posit = genomics_positions, theta.0 = theta,
  beta = beta, random.eff = error, mu.e = error_mu, sigma.ee = error_var, p0 = fp,
  p1 = tp, X = coverage, Z = covariates, binom.link = link_fct, phi = phi_v
)

test_that("the output corresponds to binomRegMethModelSim method description", {
  expect_false(length(out_binomRegMethModelSim) == 3)
  expect_equal(ncol(out_binomRegMethModelSim$S), n_positions)
  expect_equal(nrow(out_binomRegMethModelSim$S), n_samples)
  expect_equal(ncol(out_binomRegMethModelSim$Y), n_positions)
  expect_equal(nrow(out_binomRegMethModelSim$Y), n_samples)
  expect_equal(ncol(out_binomRegMethModelSim$theta), n_positions)
  expect_equal(nrow(out_binomRegMethModelSim$theta), n_samples)
})
