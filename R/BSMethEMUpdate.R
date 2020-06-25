#' @title One-step update inside the EM algorithm for fitting the functional parameter theta.0 and beta.
#'
#' @description Since the presence of conversion errors, the true methylated counts at a CpG sites are unknown. The first step inside this function is the E-step, calculating the conditional expectation of methylated counts given the observed counts. Then the second step involves maximizing the penalized likelihood using the avaialbe function \code{gam} in Package \code{mgcv}. This function fits the model where there is no random effect for individuals. Thus, the log-likelihood are sum of two parts; one is the normal approximated log-likelihood for binomial outcomes and the second part quantifies the penalization of the wiggliness of two functional parameters theta.0 and beta.
#' @param data a matrix with  \code{n*p} rows and the columns should include Y, X, Posit, Ctype and ID
#' @param pi.ij  an initial value (or value from last step) for the fitted probability \code{pi_ij} for each CpG site of each individual
#' @param p0 the probability of observing a methylated read when the underlying true status is unmethylated.
#' @param p1 the probability of observing a methylated read when the underlying true status is methylated.
#' @param n.k a vector of basis dimensions for the intercept and individual covariates. \code{n.k} specifies an upper limit of the degrees of each functional parameters.
#' @param binom.link the link function used in the binomial regression model; the default is the logit link
#' @param method the method used to estimate the smoothing parameters. The default is the "REML" method which is generally better than prediction based criterion \code{GCV.cp}
#' @param Z the covariate matrix used in BSMethEM
#' @param my.covar.fm the formula fitted in the GAM
#' @return The function returns a list of following objects
#' \itemize{
#' \item \code{pi.ij}  fitted methylation proportions
#' \item \code{par} updated basis coefficients
#' \item \code{lambda} updated smoothing parameters
#' \item \code{edf1} Effective degree of freedoms for each smooth terms in the model}
#' @author  Kaiqiong Zhao
#' @examples
#' #------------------------------------------------------------#
#' @importFrom mgcv gam
#' @export

BSMethEMUpdate <- function(data, pi.ij, p0, p1, n.k, binom.link, method, Z, my.covar.fm, Quasi = T, scale) {
  if (!(nrow(data) == length(pi.ij))) message("The row of data should be compatible with the length of initial value pi.ij")
  # The E-step
  # Calculate the "posterior" probability
  eta.1 <- p1 * pi.ij / (p1 * pi.ij + p0 * (1 - pi.ij)) # posterior probability given an observed methylated rates, what is the probability that the reads are truely methylated
  eta.0 <- (1 - p1) * pi.ij / ((1 - p1) * pi.ij + (1 - p0) * (1 - pi.ij))

  Y <- data$Y
  X <- data$X
  E.S <- Y * eta.1 + (X - Y) * eta.0

  if (Quasi) {
    gam.int.see <- suppressWarnings(mgcv::gam(as.formula(paste0("E.S/X ~", my.covar.fm)),
      family = quasibinomial(link = binom.link), weights = X,
      data = data, method = method, scale = scale
    ))
  } else {
    gam.int.see <- suppressWarnings(mgcv::gam(as.formula(paste0("E.S/X ~", my.covar.fm)),
      family = binomial(link = binom.link), weights = X,
      data = data, method = method
    ))
  }

  p_res <- residuals(gam.int.see, type = "pearson")
  d_res <- residuals(gam.int.see, type = "deviance")

  phi_fletcher <- summary(gam.int.see)$dispersion # this is actually the fixed scale paramters in the input
  out <- list(
    pi.ij = gam.int.see$fitted.values, par = gam.int.see$coefficients,
    lambda = gam.int.see$sp, edf1 = gam.int.see$edf1, pearson_res = p_res, deviance_res = d_res,
    edf = gam.int.see$edf, phi_fletcher = phi_fletcher, GamObj = gam.int.see, E.S = E.S
  )
  return(out)
}
