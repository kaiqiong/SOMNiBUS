#' @title One-step update inside the EM algorithm for fitting the
#' functional parameter theta.0 and beta.
#'
#' @description Since the presence of conversion errors, the true
#' methylated counts at a CpG sites are unknown. The first step inside this 
#' function is the E-step, calculating the conditional expectation of methylated
#' counts given the observed counts. Then the second step involves maximizing
#' the penalized likelihood using the available function.
#' \code{gam} in Package \code{mgcv}. This function fits the model where there
#' is no random effect for individuals. Thus, the log-likelihood is sum of two
#' parts; one is the normal approximated log-likelihood for binomial outcomes
#' and the second part quantifies the penalization of the wiggliness of two
#' functional parameters theta.0 and beta.
#' @param data a matrix with  \code{n*p} rows and the columns should include Y,
#' X, Posit, Ctype and ID.
#' @param pi.ij  an initial value (or value from last step) for the fitted
#' probability \code{pi_ij} for each CpG site of each individual.
#' @param p0 the probability of observing a methylated read when the underlying
#' true status is unmethylated. \code{p0} is the rate of false methylation
#' calls, i.e. false positive rate.
#' @param p1 the probability of observing a methylated read when the underlying
#' true status is methylated. \code{1-p1} is the rate of false non-methylation
#' calls, i.e. false negative rate.
#' @param n.k a vector of basis dimensions for the intercept and individual
#' covariates.
#' \code{n.k} specifies an upper limit of the degrees of each functional
#' parameters.
#' @param binom.link the link function used in the binomial regression model;
#' the default is the logit link.
#' @param method the method used to estimate the smoothing parameters. The 
#' default is the 'REML' method which is generally better than prediction based
#' criterion \code{GCV.cp}.
#' @param Z the covariate matrix used in binomRegMethModel.
#' @param my.covar.fm the formula fitted in the GAM.
#' @param Quasi whether a Quasi-likelihood estimation approach will be used; in
#' other words, whether a multiplicative dispersion is added in the model or not
#' @param scale negative values mean scale parameter should be estimated; if a
#' positive value is provided, a fixed scale will be used.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return The function returns a list of following objects
#' \itemize{
#' \item \code{pi.ij} fitted methylation proportions;
#' \item \code{par} updated basis coefficients;
#' \item \code{lambda} updated smoothing parameters;
#' \item \code{edf1} effective degrees of freedom for each smooth term in the
#' model; tr(2A - A^2), where A is the pseudo-hat matrix;
#' \item \code{pearson_res} Pearson's residuals;
#' \item \code{deviance_res} deviance residuals;
#' \item \code{edf} effective degrees of freedom for each smooth term in the
#' model; tr(A);
#' \item \code{phi_fletcher} Fletcher-based dispersion estimate;
#' \item \code{GamObj} an updated Gam output;
#' \item \code{E.S} expected methylated counts calculated from the E step.
#' }
#' @author  Kaiqiong Zhao
#' @importFrom mgcv gam
#' @importFrom stats quasibinomial residuals binomial
#' 
#' @noRd
binomRegMethModelUpdate <- function(data, pi.ij, p0, p1, n.k, binom.link,
                                    method, Z, my.covar.fm, Quasi = TRUE, 
                                    scale, reml.scale, verbose = TRUE) {
  
  t0 <- Sys.time()
  msg <- paste("One-step update inside the EM algorithm",
               "for fitting the functional parameters",
               "theta.0 and beta.")
  if(verbose) Message(msg)
  
  binomRegMethModelUpdateChecks(data = data, pi.ij = pi.ij)
  
  ## The E-step Calculate the 'posterior' probability posterior
  ## probability given an observed methylated rates, what is the
  ## probability that the reads are truely methylated
  eta.1 <- p1 * pi.ij/(p1 * pi.ij + p0 * (1 - pi.ij))
  eta.0 <- (1 - p1) * pi.ij/((1 - p1) * pi.ij + (1 - p0) * (1 - pi.ij))
  
  Y <- data$Y
  X <- data$X
  E.S <- Y * eta.1 + (X - Y) * eta.0
  
  if (Quasi) {
    gam.int.see <- suppressWarnings(mgcv::gam(as.formula(paste0("E.S/X ~",
                                                                my.covar.fm)), 
                                              family = quasibinomial(link = binom.link), 
                                              weights = X,
                                              data = data, 
                                              method = method, scale = scale))
  } else {
    gam.int.see <- suppressWarnings(mgcv::gam(as.formula(paste0("E.S/X ~",
                                                                my.covar.fm)), 
                                              family = binomial(link = binom.link), 
                                              weights = X,
                                              data = data, method = method))
  }
  
  p_res <- residuals(gam.int.see, type = "pearson")
  d_res <- residuals(gam.int.see, type = "deviance")
  ## this is actually the fixed scale paramters in the input
  data$Y <- E.S
  phi_fletcher <- phiFletcher(data, Quasi, reml.scale = reml.scale, 
                              scale = scale,
                              gam.int = gam.int.see, verbose = FALSE)
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(out <- list(pi.ij = gam.int.see$fitted.values, 
                     par = gam.int.see$coefficients,
                     lambda = gam.int.see$sp, 
                     edf1 = gam.int.see$edf1, 
                     pearson_res = p_res,
                     deviance_res = d_res, 
                     edf = gam.int.see$edf, 
                     phi_fletcher = phi_fletcher,
                     GamObj = gam.int.see, 
                     E.S = E.S))
}

#' @title Some checks for binomRegMethModelUpdate
#'
#' @description Check if inputs fit one anoter according
#' to there shapes
#' @param data a matrix with  \code{n*p} rows and the columns
#' should include Y,
#' X, Posit, Ctype and ID
#' @param pi.ij  an initial value (or value from last step) for
#' the fitted probability
#' \code{pi_ij} for each CpG site of each individual
#' @author SLL
#' @noRd
binomRegMethModelUpdateChecks <- function(data, pi.ij) {
  if (!(nrow(data) == length(pi.ij))) {
    e_msg <- paste("The row of data should be compatible with the length",
                   "of initial value pi.ij")
    Error(e_msg)
  }
}