#' @title A smoothed-EM algorithm to estimate covariate effects and test regional association in Bisulfite Sequencing-derived methylation data
#'
#' @description This function fits a binomial regression model where the outcome - methylated reads- are contaminated by known error rates \code{p0} and \code{p1} and the covariate effects are smoothly varying across genomic positions. The functional parameters of the smooth covariate effects are first represented by a linear combination of a bunch of restricted cubic splines (with dimention \code{n.k}), and a smoothness penalization term which depends on the smoothing parameters \code{lambdas} is also added to control smoothness.
#' @description The estimation is performed by an iterated EM algorithm. Each M step constitutes an outer Newton's iteration to estimate smoothing parameters \code{lambdas} and an inner P-IRLS iteration to estimate spline coefficients \code{alpha} for the covariate effects. Currently, the computation in the M step depends the implementation of \code{gam()} in package \code{mgcv}.
#' @param data a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site) and `ID` (sample ID). The covariate information, such as disease status or cell type composition are listed in column 5 and onwards.
#' @param n.k a vector of basis dimensions for the intercept and individual covariates. \code{n.k} specifies an upper limit of the degrees of each functional parameters.
#' @param p0 the probability of observing a methylated read when the underlying true status is unmethylated.
#' @param p1 the probability of observing a methylated read when the underlying true status is methylated.
#' @param covs a vector of covariate names. The covariates with names in \code{covs} will be included in the model and their covariate effects will be estimated.
#' @param Quasi whether a Quasi-likelihood estimation approach will be used
#' @param epsilon numeric; stopping criterion for the closeness of estimates of spline coefficients from two consecutive iterations.
#' @param epsilon.lambda numeric; stopping criterion for the closeness of estimates of smoothing parameter \code{lambda} from two consecutive iterations.
#' @param maxStep the algorithm will step if the iteration steps exceed \code{maxStep}
#' @param detail indicate whether print the number of iterations
#' @param binom.link the link function used in the binomial regression model; the default is the logit link
#' @param method the method used to estimate the smoothing parameters. The default is the "REML" method which is generally better than prediction based criterion \code{GCV.cp}
#' @param RanEff whether sample-level random effects are added or not
#' @param reml.scale whether a REML-based scale (dispersion) estimator is used. The default is Fletcher-based estimator
#' @param scale nagative values mean scale paramter should be estimated; if a positive value is provided, a fixed scale will be used.
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{est}: estimates of the spline basis coefficients alphas
#' \item \code{lambda}: estimates of the smoothing parameters for each functional paramters
#' \item \code{est.pi}: predicted methylation lelves for each row in the input \code{data}
#' \item \code{ite.points}: estimates of \code{est}, \code{lambda} at each EM iteration
#' \item \code{cov1}: estimated variance-covariance matrix of the basis coefficients alphs
#' \item \code{reg.out}: output from the reginal zero effect test
#' \item \code{chi.sq}: chi-square statistics for each covariates (including intercept) w.r.t the zero regional effect tests
#' \item \code{pvalue}: the p value from the regional tests for each covariate
#' \item \code{pvalue.log}: the log of \code{pvalue}
#' \item \code{SE.out} a matrix of the estimated pointwise Standard Errors (SE); number of rows are the number of unique CpG sites in the input data and the number of columns equal to the total number of covariates fitted in the model (the first one is the intercept)
#' \item \code{SE.pos} the genomic postions for each row of CpG sites in the matrix \code{SE.out}
#' \item \code{Beta.out} a matrix of the estimated covariate effects beta(t), here t denots the genomic positions.
#' \item \code{ncovs} number of functional paramters in the model (including the intercept)
#' }
#' @author  Kaiqiong Zhao
#' @seealso  \link[mgcv]{gam}
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' out <- BSMethEM(
#'   data = RAdat.f, n.k = rep(5, 3), p0 = 0.003307034, p1 = 0.9,
#'   epsilon = 10^(-6), epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE
#' )
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom mgcv model.matrix.gam
#' @importFrom mgcv s
#' @importFrom Matrix bdiag
#' @importFrom stats as.formula binomial pchisq rbinom rnorm
#' @export
BSMethEM <- function(data, n.k, p0 = 0.003, p1 = 0.9, Quasi = TRUE, epsilon = 10^(-6), epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE, binom.link = "logit", method = "REML", covs = NULL, RanEff = T, reml.scale = F, scale = -2) {
  n.k <<- n.k # an error of 'n.k is not found' would appear if without this golable environment assignment; so I save n.k in a parent scope
  data <- data.frame(data)

  if (is.factor(data$Position)) {
    # message("The Position in the data set should be numeric other than a factor")
    data$Position <- as.numeric(as.character(data$Position))
  }

  if (any(!c("Meth_Counts", "Total_Counts", "Position") %in% colnames(data))) stop('Please make sure object "data" have columns named as "Meth_Counts", "Total_Counts" and "Position" ')


  colnames(data)[match(c("Meth_Counts", "Total_Counts", "Position"), colnames(data))] <- c("Y", "X", "Posit")

  if (is.null(covs)) {
    Z <- as.matrix(data[, -(1:4)], ncol = ncol(data) - 4)
    colnames(Z) <- colnames(data)[-c(1:4)]
  } else {
    id <- match(covs, colnames(data))
    if (any(is.na(id))) {
      stop(paste(covs[is.na(id)], " is not found in the input data frame"))
    } else {
      Z <- as.matrix(data[, id], ncol = length(id))
      colnames(Z) <- covs
    }
  }

  if (length(n.k) != (ncol(Z) + 1)) stop("The length of n.k should equal to the number of covariates plus 1 (for the intercept)")
  if (any(data$X == 0)) stop("The rows with Total_Counts equal to 0 should be deleted beforehand")
  if (any(is.na(Z))) stop("The covariate information should not have missing values")
  if (any(!is.numeric(Z))) stop("Please transform the covariate information into numeric values, eg. use dummy variables for the categorical covariates")

  # The smoothing formula corresponding to the Z
  formula.z.part <- sapply(1:(ncol(Z)), function(i) {
    paste0("s(Posit, k = n.k[", i + 1, "], fx=F, bs=\"cr\", by = Z[,", i, "])")
  })
  my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=F, bs=\"cr\")", formula.z.part), collapse = "+")
  if (RanEff) {
    my.covar.fm <- paste0(my.covar.fm, "+ s(ID, bs = \"re\")")
    data$ID <- as.factor(data$ID)
  }
  # Fit gam for the initial value
  if (Quasi) {
    gam.int <- mgcv::gam(as.formula(paste0("Y/X ~", my.covar.fm)), family = quasibinomial(link = binom.link), weights = X, data = data, method = method, scale = scale)
  } else {
    gam.int <- mgcv::gam(as.formula(paste0("Y/X ~", my.covar.fm)), family = binomial(link = binom.link), weights = X, data = data, method = method, scale = scale)
  }
  # Estimates
  old.pi.ij <- gam.int$fitted.values
  old.par <- gam.int$coefficients
  lambda <- gam.int$sp
  # phi_fletcher = summary(gam.int)$dispersion
  p_res <- residuals(gam.int, type = "pearson")
  d_res <- residuals(gam.int, type = "deviance")
  edf.out <- gam.int$edf
  edf1.out <- gam.int$edf1

  # Note: this phi_fletcher can be also self-calculated
  if (Quasi & scale <= 0) { # calculate the estimate of phi if Quasi = T and scale is unknown
    my_s <- (1 - 2 * old.pi.ij) / (data$X * old.pi.ij * (1 - old.pi.ij)) * (data$Y - data$X * old.pi.ij) #* sqrt(data$X)

    # Note: the estimator implemented in the mgcv calculated my_s with an additional multiplier sqrt(data$X)
    # But from the paper there shouldn't be this one
    # phi_p = sum( p_res^2)/(length(data$Y) - sum(edf.out))

    # Feb 21, 2020  use edf1 to calculate the pearson's dispersion estimate not the edf; because I will use the asympototic chi-square dist. of phi.est
    # March 3, 2020, the dispersion parameter estimated in gam use the edf.out instead of edf1.out
    phi_p <- sum(p_res^2) / (length(data$Y) - sum(edf.out))

    phi_fletcher <- phi_p / (1 + mean(my_s))
    if (reml.scale) {
      phi_fletcher <- gam.int$reml.scale
    }
  }

  if (!Quasi) {
    phi_fletcher <- 1
  }

  if (scale > 0) {
    phi_fletcher <- scale
  }

  if (p0 == 0 & p1 == 1) {
    out <- list(
      pi.ij = gam.int$fitted.values, par = gam.int$coefficients,
      lambda = gam.int$sp, edf1 = gam.int$edf1, pearson_res = p_res, deviance_res = d_res,
      edf = gam.int$edf, phi_fletcher = phi_fletcher, GamObj = gam.int
    )
    new.par <- out$par
    new.lambda <- out$lambda
    new.pi.ij <- out$pi.ij

    Est.points <- c(new.par, new.lambda, phi_fletcher)
  } else {
    out <- BSMethEMUpdate(data, old.pi.ij, p0 = p0, p1 = p1, n.k = n.k, binom.link = binom.link, method = method, Z = Z, my.covar.fm = my.covar.fm, Quasi = Quasi, scale = phi_fletcher)
    new.par <- out$par
    new.lambda <- out$lambda
    new.pi.ij <- out$pi.ij
    new.phi <- out$phi_fletcher
    i <- 1
    Est.points <- rbind(c(old.par, lambda, phi_fletcher), c(new.par, new.lambda, new.phi))
    # Do the iteration
    # The stopping criterion is that estimator a and b are close enough
    # I exclude the criterio that lambda-a and lambda-b are close

    while (sqrt(sum((new.pi.ij - old.pi.ij)^2)) > epsilon & i < maxStep) {
      i <- i + 1
      old.par <- new.par
      old.pi.ij <- new.pi.ij

      out <- BSMethEMUpdate(data, old.pi.ij, p0 = p0, p1 = p1, binom.link = binom.link, method = method, Z = Z, my.covar.fm = my.covar.fm, Quasi = Quasi, scale = phi_fletcher)
      new.par <- out$par
      new.lambda <- out$lambda
      new.pi.ij <- out$pi.ij
      new.phi <- out$phi_fletcher

      Est.points <- rbind(Est.points, c(new.par, new.lambda, new.phi))
      if (detail) {
        print(paste0("iteration", i))
      }
    }
    # Update
    phi_fletcher <- out$phi_fletcher
  }
  # Effective degrees of freedom:  edf1 -- good for chisquare test and p value calculation tr(2A - A^2)
  edf1.out <- out$edf1
  # Effective degree of freedom: edf --trace of the hat matrix
  edf.out <- out$edf
  # the residuals degrees of freedom
  resi_df <- nrow(data) - sum(edf.out)
  # Pearson Residuals
  p_res <- out$pearson_res
  # Estimated dispersion paramters (Fletcher adjusted)
  phi_fletcher <- out$phi_fletcher

  GamObj <- out$GamObj



  #--------------------------------------------
  # Calculate SE of alpha
  #-------------------------------------------
  # the part to calculate the SE of alpha
  my.design.matrix <- mgcv::model.matrix.gam(GamObj) # the model matrix for the GamObj and the FinalGamObj is the same. the difference is only on the outcomes
  Y <- data$Y
  X <- data$X
  pred.pi <- new.pi.ij
  N <- length(unique(data$ID))
  #---------------- this part is inside  the trycatch

  phi_reml <- GamObj$reml.scale
  # phi_gam_default <- GamObj$scale
  #-------------------------------------------------------------------------
  # from the variance-covariance of alpha, to report
  # 1. var(beta_p(t))
  # 2. Report the chi-square statistic for each covariate
  # 3. Report the p-value for each covariate
  #-------------------------------------------------------------------------

  # ------ 3: estimate of beta(t) --- #
  uni.pos <- unique(data$Posit)
  uni.id <- match(uni.pos, data$Posit)
  BZ <- my.design.matrix[uni.id, 1:n.k[1]]
  BZ.beta <- lapply(1:ncol(Z), function(i) {
    mgcv::smooth.construct(mgcv::s(Posit,
      k = n.k[i + 1],
      fx = F, bs = "cr"
    ),
    data = data[uni.id, ],
    knots = NULL
    )$X
  })

  # Use PredictMat to get BZ.beta

  cum_s <- cumsum(n.k)
  alpha.sep <- lapply(1:ncol(Z), function(i) {
    new.par[(cum_s[i] + 1):cum_s[i + 1]]
  })
  alpha.0 <- new.par[1:n.k[1]]

  Beta.out <- cbind(BZ %*% alpha.0, sapply(1:ncol(Z), function(i) {
    BZ.beta[[i]] %*% alpha.sep[[i]]
  }))

  colnames(Beta.out) <- c("Intercept", colnames(Z))

  # Simply check to see
  # tryCatch({
  #----------------------------------------------------------------
  # calculate var_cov (for alpha & beta)  from the Hessian matrix
  #----------------------------------------------------------------

  H <- Hessian(w_ij = pred.pi * (1 - pred.pi) / phi_fletcher, new.par, new.lambda, X, Y, my.design.matrix, GamObj, Z, pred.pi, p0, p1, disp_est = phi_fletcher, RanEff = RanEff, N = N)

  var.cov.alpha <- solve(-H)
  var.alpha.0 <- var.cov.alpha[1:n.k[1], 1:n.k[1]]
  var.alpha.sep <- lapply(1:ncol(Z), function(i) {
    var.cov.alpha[(cum_s[i] + 1):cum_s[i + 1], (cum_s[i] + 1):cum_s[i + 1]]
  })
  # solve and Ginv give very similar results, but MASS is very different from the other two



  # A more efficient way to calculate SE rowsum(A*B) is faster than diag(A %*% t(B))
  SE.out <- cbind(
    sqrt(pmax(0, rowSums((BZ %*% var.alpha.0) * BZ))),
    sapply(1:ncol(Z), function(i) {
      sqrt(pmax(0, rowSums((BZ.beta[[i]] %*% var.alpha.sep[[i]]) * BZ.beta[[i]])))
    })
  )

  rownames(SE.out) <- uni.pos
  colnames(SE.out) <- c("Intercept", colnames(Z))
  SE.pos <- uni.pos



  SE.out.REML.scael <- SE.out / sqrt(phi_fletcher) * sqrt(phi_reml) # phi_reml

  #----------------------------------------------------------------
  # calculate the region-based statistic from the testStats function in mgcv
  #----------------------------------------------------------------

  # A more efficient way to extract design matrix. use a random sample of rows of the data to reduce the computational cost
  if (RanEff) {
    re.test <- T
  } else {
    re.test <- F
  }

  if (!is.null(GamObj$R)) {
    X_d <- GamObj$R
  } else {
    sub.samp <- max(1000, 2 * length(GamObj$coefficients))
    if (nrow(GamObj$model) > sub.samp) { ## subsample to get X for p-values calc.
      seed <- try(get(".Random.seed", envir = .GlobalEnv), silent = TRUE) ## store RNG seed
      if (inherits(seed, "try-error")) {
        runif(1)
        seed <- get(".Random.seed", envir = .GlobalEnv)
      }
      kind <- RNGkind(NULL)
      RNGkind("default", "default")
      set.seed(11) ## ensure repeatability
      ind <- sample(1:nrow(GamObj$model), sub.samp, replace = FALSE) ## sample these rows from X
      X_d <- predict(GamObj, GamObj$model[ind, ], type = "lpmatrix")
      RNGkind(kind[1], kind[2])
      assign(".Random.seed", seed, envir = .GlobalEnv) ## RNG behaves as if it had not been used
    } else { ## don't need to subsample
      X_d <- model.matrix(GamObj)
    }
    X_d <- X_d [!is.na(rowSums(X_d)), ] ## exclude NA's (possible under na.exclude)
  } ## end if (m>0)

  # p : estimated paramters --- alpha.0, alpha.seq
  # Xt: the design matrix --- my.design.matrix
  # V: estimated variance  matrix
  ## on entry `rank' should be an edf estimate
  ## 0. Default using the fractionally truncated pinv.
  ## 1. Round down to k if k<= rank < k+0.05, otherwise up.
  ## res.df is residual dof used to estimate scale. <=0 implies
  ## fixed scale.

  s.table <- BSMethEM_summary(GamObj, var.cov.alpha, new.par, edf.out, edf1.out, X_d, resi_df, Quasi, scale, RanEff, re.test, Z)
  s.table.REML.scale <- BSMethEM_summary(GamObj, var.cov.alpha / phi_fletcher * phi_reml, new.par, edf.out, edf1.out, X_d, resi_df, Quasi, scale, RanEff, re.test, Z)
  # var_out = list(cov1 = var.cov.alpha, reg.out = reg.out,  SE.out = SE.out, uni.pos = SE.pos,  pvalue = pvalue , ncovs = ncol(Z)+1)
  # Est_out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points,
  #               Beta.out = Beta.out, phi_fletcher = phi_fletcher)


  if (RanEff) {
    sigma00 <- GamObj$reml.scale / GamObj$sp["s(ID)"]
  } else {
    sigma00 <- NA
  }

  reg.out.gam <- summary(GamObj)$s.table


  return(out = list(
    est = new.par, lambda = new.lambda, est.pi = new.pi.ij,
    Beta.out = Beta.out,
    phi_fletcher = phi_fletcher,
    phi_reml = phi_reml,
    reg.out = s.table,
    reg.out.reml.scale = s.table.REML.scale,
    cov1 = var.cov.alpha,
    reg.out.gam = reg.out.gam,
    SE.out = SE.out,
    SE.out.REML.scale = SE.out.REML.scael,
    uni.pos = SE.pos,
    ncovs = ncol(Z) + 1, ite.points = Est.points, sigma00 = sigma00
  ))
}
#' @title Hessian Lorem ipsum dolor sit amet
#'
#' @description Lorem ipsum dolor sit amet
#' @description Lorem ipsum dolor sit amet
#' @param w_ij Lorem ipsum dolor sit amet
#' @param new.par Lorem ipsum dolor sit amet
#' @param new.lambda Lorem ipsum dolor sit amet
#' @param X Lorem ipsum dolor sit amet
#' @param Y Lorem ipsum dolor sit amet
#' @param my.design.matrix Lorem ipsum dolor sit amet
#' @param gam.int Lorem ipsum dolor sit amet
#' @param Z Lorem ipsum dolor sit amet
#' @param pred.pi Lorem ipsum dolor sit amet
#' @param p0 Lorem ipsum dolor sit amet
#' @param p1 Lorem ipsum dolor sit amet
#' @param disp_est Lorem ipsum dolor sit amet
#' @param RanEff Lorem ipsum dolor sit amet
#' @param N Lorem ipsum dolor sit amet
#' @return Lorem ipsum dolor sit amet
#' \itemize{
#' \item Lorem ipsum dolor sit amet
#' }
#' @author  Kaiqiong Zhao
#' @examples
#' #------------------------------------------------------------#
#' # Lorem ipsum dolor sit amet
#' @importFrom Matrix bdiag
Hessian <- function(w_ij, new.par, new.lambda, X, Y, my.design.matrix, gam.int, Z, pred.pi, p0, p1, disp_est, RanEff, N) {

  # Q1: the second partial derivative w.r.t alpha^2
  # Q2: the second derivative w.r.t alpha & alpha_star
  res <- outer(1:length(new.par), 1:length(new.par), Vectorize(function(l, m) sum(-X * w_ij * my.design.matrix[, m] * my.design.matrix[, l])))
  smoth.mat <- lapply(as.list(1:(ncol(Z) + 1)), function(i) {
    gam.int$smooth[[i]]$S[[1]] * new.lambda[i]
  }) # extract the penalty matrix

  smoth.mat[[length(smoth.mat) + 1]] <- 0 # assume the lambda for the constant of the intercept is 0 -- no penalization
  if (RanEff) {
    smoth.mat[[length(smoth.mat) + 1]] <- diag(N) * new.lambda[ncol(Z) + 2] # !!!! Otherwise, we get very wide CI
    span.penal.matrix <- as.matrix(Matrix::bdiag(smoth.mat[c(length(smoth.mat) - 1, (1:(length(smoth.mat) - 2)), length(smoth.mat))]))
  } else {
    span.penal.matrix <- as.matrix(Matrix::bdiag(smoth.mat[c(length(smoth.mat), (1:(length(smoth.mat) - 1)))]))
  }

  Q1_with_lambda <- res - span.penal.matrix / disp_est
  Q1_no_lambda <- res

  Q2 <- outer(1:length(new.par), 1:length(new.par), Vectorize(function(l, m) {
    term1 <- Y * p1 * p0 / (p1 * pred.pi + p0 * (1 - pred.pi))^2 + (X - Y) * (1 - p1) * (1 - p0) / ((1 - p1) * pred.pi + (1 - p0) * (1 - pred.pi))^2
    sum(term1 * w_ij * my.design.matrix[, m] * my.design.matrix[, l])
  }))

  return(Q1_with_lambda + Q2)
}
#' @title BSMethEM_summary Lorem ipsum dolor sit amet
#'
#' @description Lorem ipsum dolor sit amet
#' @description Lorem ipsum dolor sit amet
#' @param GamObj Lorem ipsum dolor sit amet
#' @param var.cov.alpha Lorem ipsum dolor sit amet
#' @param new.par Lorem ipsum dolor sit amet
#' @param edf.out Lorem ipsum dolor sit amet
#' @param edf1.out Lorem ipsum dolor sit amet
#' @param X_d Lorem ipsum dolor sit amet
#' @param gam.int Lorem ipsum dolor sit amet
#' @param resi_df Lorem ipsum dolor sit amet
#' @param Quasi Lorem ipsum dolor sit amet
#' @param scale Lorem ipsum dolor sit amet
#' @param p1 Lorem ipsum dolor sit amet
#' @param disp_est Lorem ipsum dolor sit amet
#' @param RanEff Lorem ipsum dolor sit amet
#' @param re.test Lorem ipsum dolor sit amet
#' @param Z Lorem ipsum dolor sit amet
#' @return Lorem ipsum dolor sit amet
#' \itemize{
#' \item Lorem ipsum dolor sit amet
#' }
#' @author  Kaiqiong Zhao
#' @examples
#' #------------------------------------------------------------#
#' # Lorem ipsum dolor sit amet
#' @import mgcv
BSMethEM_summary <- function(GamObj, var.cov.alpha, new.par, edf.out, edf1.out, X_d, resi_df, Quasi, scale, RanEff, re.test, Z) {
  ii <- 0
  m <- length(GamObj$smooth)
  df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
  for (i in 1:m) {
    start <- GamObj$smooth[[i]]$first.para
    stop <- GamObj$smooth[[i]]$last.para

    V <- var.cov.alpha[start:stop, start:stop, drop = FALSE] ## Bayesian

    p <- new.par[start:stop] # params for smooth
    edfi <- sum(edf.out[start:stop]) # edf for this smooth
    ## extract alternative edf estimate for this smooth, if possible...
    edf1i <- sum(edf1.out[start:stop])
    Xt <- X_d [, start:stop, drop = FALSE]
    fx <- if (inherits(GamObj$smooth[[i]], "tensor.smooth") &&
      !is.null(GamObj$smooth[[i]]$fx)) {
      all(GamObj$smooth[[i]]$fx)
    } else {
      GamObj$smooth[[i]]$fixed
    }
    if (!fx && GamObj$smooth[[i]]$null.space.dim == 0 && !is.null(GamObj$R)) { ## random effect or fully penalized term
      res <- if (re.test) mgcv:::reTest(GamObj, i) else NULL # Test the mth smooth for equality to zero  (m is not the RE term)
      ## and accounting for all random effects in model
    } else { ## Inverted Nychka interval statistics

      if (Quasi) rdf <- resi_df else rdf <- -1
      res <- mgcv:::testStat(p, Xt, V, min(ncol(Xt), edf1i), type = 0, res.df = rdf)
    }

    if (!is.null(res)) {
      ii <- ii + 1
      df[ii] <- res$rank
      chi.sq[ii] <- res$stat
      s.pv[ii] <- res$pval
      edf1[ii] <- edf1i
      edf[ii] <- edfi
      names(chi.sq)[ii] <- GamObj$smooth[[i]]$label
    }

    if (ii == 0) {
      df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, 0)
    } else {
      df <- df[1:ii]
      chi.sq <- chi.sq[1:ii]
      edf1 <- edf1[1:ii]
      edf <- edf[1:ii]
      s.pv <- s.pv[1:ii]
    }
    if (!Quasi || scale >= 0) {
      s.table <- cbind(edf, df, chi.sq, s.pv)
      dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
    } else {
      s.table <- cbind(edf, df, chi.sq / df, s.pv)
      dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "F", "p-value"))
    }
  }
  s.table <- s.table[, -1]
  if (RanEff) {
    rownames(s.table) <- c("Intercept", colnames(Z), "ID")
  } else {
    rownames(s.table) <- c("Intercept", colnames(Z))
  }
  colnames(s.table)[1] <- "EDF"
  s.table
}
