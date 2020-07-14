#' @title Simulate Bisulfite sequencing data from specified smooth covariate effects
#'
#' @description Simulate Bisulfite sequencing data from a Generalized Additive Model with functional parameters varying with the genomic position. Both the true methylated counts and observed methylated counts are generated, given the error/conversion rate parameters \code{p0} and \code{p1}.
#'
#' @param n sample size
#' @param posit  genomic position; a numeric vector of size \code{p} (the number of CpG sites in the considered region).
#' @param theta.0 a functional parameter for the intercept of the GAMM model; a numeric vector of size \code{p}.
#' @param random.eff indicate whether adding the subject-specific random effect term \code{e}.
#' @param mu.e the mean of the random effect; a single number.
#' @param sigma.ee variance of the random effect; a single positive number.
#' @param beta a functional parameter for the slope of cell type composition. a numeric vector of size \code{p}
#' @param p0 the probability of observing a methylated read when the underlying true status is unmethylated.
#' @param p1 the probability of observing a methylated read when the underlying true status is methylated.
#' @param X the matrix of the read coverage for each CpG in each sample; a matrix of n rows and \code{p} columns
#' @param Z numeric vector of length p for the covariate; currently, the covariate is the percentage of Cell type A (considering that the samples are composed of two cell types A and B); (Added on Feb 2018), Z can be a matrix; we allow for more than one covariates
#' @param binom.link the link  function used for simulation
#' @return The function returns a list of following objects
#' @return \code{S} the true methylation counts; a numeric matrix of \code{n} rows and \code{p} columns
#' @return \code{Y} the observed methylation counts; a numeric matrix of \code{n} rows and \code{p} columns
#' @return \code{theta} the methylation parameter (after the logit transformation); a numeric matrix of \code{n} rows and \code{p} columns
#' @author  Kaiqiong Zhao
#' @export
BSMethSim <- function(n, posit, theta.0, beta, random.eff = FALSE, mu.e = 0, sigma.ee = 1,
                      p0 = 0.003, p1 = 0.9, X, Z, binom.link = "logit") {
  if (!is.matrix((Z))) {
    message("covariate Z is not a matrix")
  }
  # if( !is.matrix(beta) ) message ('the covariate effect parameter beta is not a
  # matrix')

  if (!(nrow(X) == nrow(Z) & nrow(X) == n)) {
    message("Both X and Z should have n rows")
  }
  if (!(ncol(X) == nrow(theta.0) & ncol(X) == nrow(beta) & ncol(X) == length(posit))) {
    message("The columns of X should be the same as length of beta theta.0 and posit; They all equals to the number of CpGs")
  }
  if (ncol(beta) != ncol(Z)) {
    message("beta and Z should have the same dimentions")
  }

  # the random effect term
  if (random.eff == TRUE) {
    my.e <- rnorm(n, mean = mu.e, sd = sqrt(sigma.ee))
  } else {
    my.e <- rep(mu.e, n)
  }

  my.theta <- t(sapply(seq_len(n), function(i) {
    theta.0 + rowSums(sapply(seq_len(ncol(Z)), function(j) {
      Z[i, j] * beta[, j]
    })) + my.e[i]
  }))


  # Transform my.theta to my.pi for each (i, j)

  my.pi <- t(sapply(seq_len(nrow(my.theta)), function(i) {
    # exp(my.theta[i,])/(1+exp(my.theta[i,]))
    binomial(link = binom.link)$linkinv(my.theta[i, ])
  }))
  # Generate S-ij based on the my.pi and my.Z
  #---------------------------------------------------#
  my.S <- my.pi
  for (i in seq_len(nrow(my.S))) {
    for (j in seq_len(ncol(my.S))) {
      my.S[i, j] <- rbinom(1, size = X[i, j], prob = my.pi[i, j])
    }
  }
  #---------------------------------------------------#
  # Generate Y-ij based on the S-ij and the error rate (1-p1) and p0
  #---------------------------------------------------#
  my.Y <- my.S
  for (i in seq_len(nrow(my.Y))) {
    for (j in seq_len(ncol(my.Y))) {
      my.Y[i, j] <- sum(rbinom(my.S[i, j], size = 1, prob = p1)) + sum(rbinom(X[
        i,
        j
      ] - my.S[i, j], size = 1, prob = p0))
    }
  }
  out <- list(S = my.S, Y = my.Y, theta = my.theta, pi = my.pi)
}
