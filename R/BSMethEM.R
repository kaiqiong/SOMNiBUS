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
#' @examples #------------------------------------------------------------#
#' data(RAdat)
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0,])
#' out <- BSMethEM(data=RAdat.f, n.k = rep(5,3), p0 = 0.003307034, p1 = 0.9,
#' epsilon = 10^(-6), epsilon.lambda = 10^(-3), maxStep = 200, detail=FALSE)
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom mgcv model.matrix.gam
#' @importFrom mgcv s
#' @importFrom Matrix bdiag
#' @importFrom stats as.formula binomial pchisq rbinom rnorm
#' @export
BSMethEM = function (data, n.k, p0 = 0.003, p1 = 0.9, Quasi = TRUE, epsilon = 10^(-6),  epsilon.lambda = 10^(-3), maxStep = 200,  detail=FALSE, binom.link = "logit",method="REML", covs = NULL){

  n.k <<-n.k # an error of 'n.k is not found' would appear if without this golable environment assignment; so I save n.k in a parent scope
  data <- data.frame(data)

  if(is.factor(data$Position) ) stop ("The Position in the data set should be numeric other than a factor")
  colnames(data)[match(c("Meth_Counts", "Total_Counts", "Position") ,colnames(data)) ] <- c("Y", "X", "Posit")

  if(is.null(covs)){
    Z <- as.matrix(data[,-(1:4)], ncol = ncol(data)-4)
    colnames(Z) <- colnames(data)[-c(1:4)]
  }else{
   id <- match(covs, colnames(data))
   if(any(is.na(id))){
     stop(paste(covs[is.na(id)], " is not found in the input data frame"))
   }else{
     Z <- as.matrix(data[,id], ncol = length(id))
     colnames(Z) <- covs
     }
  }

  if(length(n.k)!= (ncol(data) -3)) stop('The length of n.k should equal to the number of covariates plus 1 (for the intercept)')
  if(any(data$X == 0)) stop('The rows with Total_Counts equal to 0 should be deleted beforehand')
  if(any(is.na(Z))) stop('The covariate information should not have missing values')
  if(any(!is.numeric(Z))) stop('Please transform the covariate information into numeric values, eg. use dummy variables for the categorical covariates')

  # The smoothing formula corresponding to the Z
  formula.z.part <- sapply(1:(ncol(Z)), function(i){
    paste0("s(Posit, k = n.k[",i+1,"], fx=F, bs=\"cr\", by = Z[,", i, "])"  ) } )
  my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=F, bs=\"cr\")", formula.z.part), collapse="+")

  # Fit gam for the initial value
  if(Quasi){
    gam.int <- mgcv::gam(as.formula( paste0("Y/X ~", my.covar.fm)), family =quasibinomial(link = binom.link),weights=X, data = data, method = method)
  }else{
    gam.int <- mgcv::gam(as.formula( paste0("Y/X ~", my.covar.fm)), family =binomial(link = binom.link),weights=X, data = data, method = method)
  }
  # Estimates
  old.pi.ij <- gam.int$fitted.values;old.par <-gam.int$coefficients;lambda <- gam.int$sp


  if(p0==0 & p1==1){
    p_res <- residuals(gam.int, type ='pearson')
    d_res <- residuals(gam.int, type ="deviance")
    phi_fletcher = summary(gam.int)$dispersion
    out <- list( pi.ij = gam.int$fitted.values, par = gam.int$coefficients,
                 lambda = gam.int$sp, edf1 = gam.int$edf1, pearson_res = p_res, deviance_res=d_res,
                 edf=gam.int$edf, phi_fletcher= phi_fletcher, GamObj = gam.int )
    new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij
  }else{
  # Update
  out <-  BSMethEMUpdate (data, old.pi.ij, p0 = p0, p1 = p1, n.k=n.k, binom.link = binom.link, method = method, Z = Z, my.covar.fm = my.covar.fm, Quasi = Quasi)
  new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij
  i <- 1; Est.points <- rbind(c(old.par, lambda), c(new.par, new.lambda))
  # Do the iteration
  # The stopping criterion is that estimator a and b are close enough
  # I exclude the criterio that lambda-a and lambda-b are close
  while ( sum((old.par - new.par)^2) > (epsilon/15 * length(new.par)) & i < maxStep & sum((lambda-new.lambda)^2) > epsilon.lambda){
    i <- i +1;
    old.par <- new.par
    old.pi.ij <- new.pi.ij

    out <-  BSMethEMUpdate (data, old.pi.ij, p0 = p0, p1 = p1,binom.link = binom.link, method = method, Z = Z, my.covar.fm = my.covar.fm, Quasi = Quasi)
    new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij

    Est.points <- rbind(Est.points, c(new.par,new.lambda))
    if(detail){
      print(paste0("iteration", i))
    }
  }
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

  gam.int <- GamObj <- out$GamObj
  # Note: this phi_fletcher can be also self-calculated

  #my_s <- (1-2*new.pi.ij)/(data$X*new.pi.ij*(1-new.pi.ij))*(data$Y-data$X*new.pi.ij)#*sqrt(data$X)

  # Note: the estimator implemented in the mgcv calculated my_s with an additional multiplier sqrt(data$X)
  # But from the paper there shouldn't be this one
  #phi_p = sum( p_res^2)/(length(data$Y) - sum(edf.out))
  #phi_f <- phi_p/(1+mean(my_s))

  #--------------------------------------------
  # Calculate SE of alpha
  #-------------------------------------------
  # the part to calculate the SE of alpha
  my.design.matrix <-  mgcv::model.matrix.gam(gam.int) # the model matrix for the gam.int and the FinalGamObj is the same. the difference is only on the outcomes
  Y <- data$Y ; X <- data$X
  pred.pi <- new.pi.ij
  #---------------- this part is inside  the trycatch


  #-------------------------------------------------------------------------
  # from the variance-covariance of alpha, to report
  # 1. var(beta_p(t))
  # 2. Report the chi-square statistic for each covariate
  # 3. Report the p-value for each covariate
  #-------------------------------------------------------------------------

  # ------ 3: estimate of beta(t) --- #
  uni.pos <- unique(data$Posit); uni.id <- match(uni.pos, data$Posit)
  BZ <- my.design.matrix[uni.id, 1:n.k[1]]
  BZ.beta = lapply(1:ncol(Z), function(i){mgcv::smooth.construct(mgcv::s(Posit, k = n.k[i+1],
                                                                         fx = F, bs = "cr") ,
                                                                 data = data[uni.id,],
                                                                 knots = NULL)$X })

  cum_s <-cumsum(n.k)
  alpha.sep <- lapply(1:ncol(Z), function(i){new.par[ (cum_s[i]+1):cum_s[i+1]]}); alpha.0 <- new.par[1:n.k[1]]

  Beta.out<-  cbind( BZ %*% alpha.0, sapply(1:ncol(Z), function(i){BZ.beta[[i]] %*% alpha.sep[[i]]}))

  colnames(Beta.out) <- c("Intercept", colnames(Z))

  # Simply check to see
  # tryCatch({
  #----------------------------------------------------------------
  # calculate var_cov (for alpha & beta)  from the Hessian matrix
  #----------------------------------------------------------------


  H = Hessian(w_ij =pred.pi * (1-pred.pi)/phi_fletcher,        new.par, new.lambda, X, Y, my.design.matrix, gam.int, Z,pred.pi, p0, p1, disp_est = phi_fletcher)


  var.cov.alpha <- solve(-H)

  # SE.out --- pointwise standard deviation
  SE.out <- vector(mode = "list", length = (ncol(Z)+1))
  names(SE.out) <- c("Intercept", colnames(Z))

  var.alpha.0 <- var.cov.alpha[1:n.k[1], 1:n.k[1]]
  var.beta.0 <-  BZ %*% var.alpha.0 %*% t(BZ)

  var.alpha.sep <- lapply(1:ncol(Z), function(i){var.cov.alpha[ (cum_s[i]+1):cum_s[i+1],(cum_s[i]+1):cum_s[i+1]]})  # SE of the effect of Zs [beta.1(t), beta.2(t), beta.3(t) ...]
  var.beta <- lapply(1:ncol(Z), function(i){BZ.beta[[i]] %*% var.alpha.sep[[i]] %*% t(BZ.beta[[i]])})

  SE.out <- cbind(sqrt(diag(var.beta.0)), sapply(1:ncol(Z), function(i){sqrt(diag(var.beta[[i]]))}))
  rownames(SE.out) <- uni.pos; colnames(SE.out) <-  c("Intercept", colnames(Z))
  SE.pos <- uni.pos


  # -------2: The chi-square/F statistics and P values for each covariate

  chi.sq <- c( t(alpha.0) %*% solve(var.alpha.0) %*% alpha.0, sapply(1:ncol(Z), function(i){t(alpha.sep[[i]]) %*% solve(var.alpha.sep[[i]]) %*% alpha.sep[[i]]}))
  #pvalue <- 1- pchisq(chi.sq, df = n.k)
  names(chi.sq) <- c("Intercept", colnames(Z))
  edf <- edf1.out

  df = sapply(1:(ncol(Z)+1), function(i){sum(edf[(cum_s[i]-n.k[i]+1):cum_s[i]])})
  #pvalue <- 1- pchisq(chi.sq, df = df)
  #pvalue <- 1- pchisq(chi.sq, df = n.k)


 names(chi.sq) <- c("Intercept", colnames(Z))


    if (!Quasi) {
      pvalue.log <-  pchisq(chi.sq, df = df, log.p = T, lower.tail = F)
      pvalue <- exp(pvalue.log)
      s.table <- data.frame(df, chi.sq,  formatC(pvalue, format ="e", digits = 2))
      dimnames(s.table) <- list(names(chi.sq), c("Edf", "Chi.sq", "p-value"))
    }else {
      pvalue.log <-  pf(chi.sq/df, df1=df, df2=resi_df, log.p = T, lower.tail = F)
      pvalue <- exp(pvalue.log)
      s.table <- data.frame("EDF" = round(df,4), "F" = round(chi.sq/df, 4),"p-value"= formatC(pvalue, format ="e", digits = 2))
      dimnames(s.table) <- list(names(chi.sq), c("Edf", "F", "p-value"))
    }

  #s.table <- s.table[, -1]

  #colnames(s.table)[1] <- "EDF"
  #var_out = list(cov1 = var.cov.alpha, reg.out = reg.out,  SE.out = SE.out, uni.pos = SE.pos,  pvalue = pvalue , ncovs = ncol(Z)+1)
  #Est_out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points,
  #               Beta.out = Beta.out, phi_fletcher = phi_fletcher)

 #reg.out <- data.frame( "EDF" = round(df,4), "Chi.sq" = round(chi.sq, 4),
   #                     "p_value" = formatC(pvalue, format ="e", digits = 2))

  return(out=list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij,
                  Beta.out = Beta.out, phi_fletcher = phi_fletcher,
                  reg.out = s.table,
                  cov1 = var.cov.alpha,  SE.out = SE.out, uni.pos = SE.pos,  pvalue = pvalue , ncovs = ncol(Z)+1, GamObj=GamObj ))
  # },
  #  error = function(err){
  #    return (out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points,
  #                       cov1 = diag(NA, length(new.par)), Beta.out = Beta.out , ncovs = ncol(Z)+1) )
  # }
  # )

}
Hessian <- function(w_ij, new.par, new.lambda, X, Y, my.design.matrix, gam.int, Z,pred.pi, p0, p1, disp_est){

  # Q1: the second partial derivative w.r.t alpha^2
  # Q2: the second derivative w.r.t alpha & alpha_star
  res <- outer( 1:length(new.par), 1:length(new.par), Vectorize(function(l,m)  sum(-X * w_ij * my.design.matrix[, m]*my.design.matrix[,l] )) )
  smoth.mat <-lapply(as.list(1:(ncol(Z)+1)), function(i){gam.int$smooth[[i]]$S[[1]] * new.lambda[i]})  # extract the penalty matrix

  smoth.mat[[length(smoth.mat) + 1]] <- 0  # assume the lambda for the constant of the intercept is 0 -- no penalization
  span.penal.matrix <- as.matrix(Matrix::bdiag( smoth.mat[c(length(smoth.mat), (1:(length(smoth.mat)-1)))] ))
  Q1_with_lambda <- res - span.penal.matrix/disp_est
  Q1_no_lambda <- res

  Q2 <- outer(1:length(new.par), 1:length(new.par), Vectorize(function(l,m){
    term1 <- Y*p1*p0/(p1*pred.pi + p0 * (1-pred.pi))^2 + (X-Y)*(1-p1)*(1-p0)/((1-p1)*pred.pi + (1-p0) * (1-pred.pi))^2
    sum(term1 * w_ij * my.design.matrix[,m] * my.design.matrix[,l])
  }))

  return(Q1_with_lambda + Q2)

}
