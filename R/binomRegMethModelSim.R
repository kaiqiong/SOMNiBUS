#' @title Simulate Bisulfite sequencing data from specified smooth covariate
#' effects
#'
#' @description Simulate Bisulfite sequencing data from a Generalized Additive
#' Model with functional parameters varying with the genomic position. Both the
#' true methylated counts and observed methylated counts are generated, given
#' the error/conversion rate parameters \code{p0} and \code{p1}. In addition,
#' the true methylated counts can be simulated from a binomial or a dispersed
#' binomial distribution (Beta-binomial distribution).
#'
#' @param n sample size
#' @param posit a numeric vector of size \code{p} (the number
#' of CpG sites in the considered region) containing the genomic positions; 
#' @param theta.0 numeric vector of size \code{p} which is a functional
#' parameter for the intercept of the GAMM model.
#' @param random.eff indicates whether adding the subject-specific random effect
#' term \code{e}.
#' @param mu.e number, the mean of the random effect.
#' @param sigma.ee positive number, variance of the random effect
#' @param beta numeric vector of size \code{p} which is a functional parameter
#' for the slope of cell type composition.
#' @param phi a vector of length \code{p} determining the multiplicative 
#' dispersion parameter for each loci in a region. The dispersed-Binomial
#' counts are simulated from beta-binomial distribution, so each element of phi
#' has to be greater than 1.
#' @param p0 the probability of observing a methylated read when the underlying
#' true status is unmethylated. \code{p0} is the rate of false methylation 
#' calls, i.e. false positive rate.
#' @param p1 the probability of observing a methylated read when the underlying
#' true status is methylated. \code{1-p1} is the rate of false non-methylation
#' calls, i.e. false negative rate.
#' @param X the matrix of the read coverage for each CpG in each sample; a 
#' matrix of n rows and \code{p} columns.
#' @param Z numeric matrix with \code{p} columns and \code{n} rows storing the
#' covariate information.
#' @param binom.link the link function used for simulation
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return The function returns a list of following objects
#' \itemize{
#' \item \code{S} a numeric matrix of \code{n} rows and \code{p} columns
#' containing the true methylation counts; 
#' \item \code{Y} a numeric matrix of \code{n} rows and \code{p} columns
#' containing the observed methylation counts; 
#' \item \code{theta} a numeric matrix of \code{n} rows and \code{p} columns
#' containing the methylation parameter (after the logit transformation);
#' \item \code{pi} a numeric matrix of \code{n} rows and \code{p} columns 
#' containing the true methylation proportions used to simulate the data.
#' }
#' @author  Kaiqiong Zhao
#' @importFrom VGAM rbetabinom
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' out <- binomRegMethModel(
#'    data=RAdat.f, n.k=rep(5, 3), p0=0, p1=1,
#'    epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200, RanEff = FALSE
#' )
#' Z = as.matrix(RAdat.f[match(unique(RAdat.f$ID), RAdat.f$ID),
#' c('T_cell', 'RA')])
#' set.seed(123)
#' X = matrix(sample(80, nrow(Z)*length(out$uni.pos), replace = TRUE),
#' nrow = nrow(Z), ncol = length(out$uni.pos))+10
#' simdat = binomRegMethModelSim(n=nrow(Z), posit= out$uni.pos,
#' theta.0=out$Beta.out[,1], beta= out$Beta.out[,-1], random.eff=FALSE,
#' mu.e=0,sigma.ee=1, p0=0.003, p1=0.9,X=X , Z=Z, binom.link='logit',
#' phi = rep(1, length(out$uni.pos)))
#' @export
binomRegMethModelSim <- function(n, posit, theta.0, beta, phi, 
                                 random.eff = FALSE, mu.e = 0, sigma.ee = 1, 
                                 p0 = 0.003, p1 = 0.9, X, Z, 
                                 binom.link = "logit", verbose = TRUE) {
  
  t0 <- Sys.time()
  msg <- paste("Simulate Bisulfite sequencing data from",
               "specified smooth covariate effects")
  if(verbose) Message(msg)
  
  ## some checks on inputs
  binomRegMethModelSimChecks(n = n, posit = posit, Z = Z, X = X, 
                             theta.0 = theta.0, beta = beta, verbose = verbose)
  
  ## the random effect term
  if (random.eff == TRUE) {
    my.e <- rnorm(n, mean = mu.e, sd = sqrt(sigma.ee))
  } else {
    my.e <- rep(mu.e, n)
  }
  my.theta <- t(vapply(seq_len(n), function(i) {
    theta.0 + rowSums(vapply(seq_len(ncol(Z)), function(j) {
      Z[i, j] * beta[, j]
    }, rep(1, length(posit)))) + my.e[i]
  }, rep(1, length(posit))))
  ## Transform my.theta to my.pi for each (i, j)
  my.pi <- t(vapply(seq_len(nrow(my.theta)), function(i) {
    ## exp(my.theta[i,])/(1+exp(my.theta[i,]))
    binomial(link = binom.link)$linkinv(my.theta[i, ])
  }, rep(1, length(posit))))
  ## Generate S-ij based on the my.pi and my.Z
  my.S <- my.pi
  for (i in seq_len(nrow(my.S))) {
    for (j in seq_len(ncol(my.S))) {
      my.S[i, j] <- VGAM::rbetabinom(1, size = X[i, j], prob = my.pi[i,j],
                                     rho = (phi[j] - 1)/(X[i, j] - 1))
    }
  }
  ## Generate Y-ij based on the S-ij and the error rate (1-p1) and p0
  my.Y <- my.S
  for (i in seq_len(nrow(my.Y))) {
    for (j in seq_len(ncol(my.Y))) {
      my.Y[i, j] <- sum(rbinom(my.S[i, j], size = 1, prob = p1)) +
        sum(rbinom(X[i, j] - my.S[i, j], size = 1, prob = p0))
    }
  }
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  out <- list(S = my.S, Y = my.Y, theta = my.theta, pi = my.pi)
}

#' @title Some checks for binomRegMethModelSim
#'
#' @description Check if inputs fit one anoter according to there shapes
#' @param n sample size
#' @param posit  genomic position; a numeric vector of size
#' \code{p} (the number of
#' CpG sites in the considered region).
#' @param Z numeric matrix with \code{p} columns and \code{n}
#' rows storing the covariate
#' information
#' @param X the matrix of the read coverage for each CpG in each
#' sample; a matrix of
#' n rows and \code{p} columns
#' @param theta.0 a functional parameter for the intercept of the
#' GAMM model; a numeric
#' vector of size \code{p}.
#' @param beta a functional parameter for the slope of cell type
#' composition. a numeric
#' vector of size \code{p}
#' @author  Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
binomRegMethModelSimChecks <- function(n, posit, Z, X, theta.0, beta, 
                                       verbose = TRUE) {
  
  t0 <- Sys.time()
  if(verbose) Message("Check inputs consistency")
  
  if (!is.matrix((Z))) {
    Error("Covariate Z is not a matrix")
  }
  if (!(nrow(X) == nrow(Z) & nrow(X) == n)) {
    Error("Both X and Z should have n rows")
  }
  if (!(ncol(X) == length(theta.0) & 
        ncol(X) == nrow(beta) & 
        ncol(X) == length(posit))) {
    e_msg <- paste("The columns of X should be the same",
                   "as length of beta theta.0 and posit;",
                   "They all equals to the number of CpGs")
    Error(e_msg)
  }
  if (ncol(beta) != ncol(Z)) {
    Error("beta and Z should have the same dimentions")
  }
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
}
