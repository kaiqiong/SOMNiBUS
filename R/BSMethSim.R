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
#' @examples #------------------------------------------------------------#
#'#Generate functional parameter theta.0, position based from the
#'#real data set "Data_HHM" in the package SMSC
#'#------------------------------------------------------------#
#' library(SMSC)
#' data("Data_HHM")
#' y <- Data_HHM
#'fit1 <- SMSC(y$Ccount, y$CT, y$position, method="KNN2")
#'par(mfrow=c(1,1))
#'plotSMSC(y$Ccount/y$CT, y$origin, y$position, fit1) # The methylation levels are stored in fit$pi
#'plot(y$position, fit1$pi, type="l", lwd=3)
# #retrieve the region between (590kb - 700kb)
#'abline(v = 590000)
#'abline(v = 700000)
#'window.id <- which( y$position > 590000 & y$position <700000)
#'# smooth.spline for the logit scale
#'my.exp.y <- log(fit1$pi[window.id]/(1-fit1$pi[window.id]) )
#'my.nknots = 5
#'plot(y$position[window.id], my.exp.y, ylim = c(-4,8), main="Region of interest",
#'    xlab = "Position", ylab = expression( paste("Logit of methylation level ", theta[0], sep=" ")) )
#'ff <- smooth.spline(y$position[window.id], my.exp.y, nknots=my.nknots)
#'lines(ff, col=my.nknots, lty=6)
#'# transform back to the scale (0, 1)
#'my.fit.pi <-  exp(ff$y)/(1+exp(ff$y))
#'
#'plot(y$position[window.id], fit1$pi[window.id], ylim = c(-0.2, 1.2), main="Region of interest",
#'     xlab = "Position", ylab = expression( paste("A smooth estimate of ", pi[0], sep=" ")))
#'lines(y$position[window.id], my.fit.pi,col = 2, type="b", lty=3 )
#'#---------------------------------------------------------------------
#'# Note: the estimation of the smooth curve pi0[t] is obtained from function
#'# "smooth.spline" -- fixing the degree of freedom (i.e. fixing lambda) -- selected by eyes
#'my.t <- y$position[window.id] # Position
#'theta.0 <- predict(ff, my.t)$y # theta.0  There are 90 CpGs in this region
#'pi.0 <- exp(theta.0)/(1+exp(theta.0))
#'#---------------------------------------------------#
#'# Sample size and random effect
#'#---------------------------------------------------#
#'my.n <- 8   # sample size
#'sigma.ee <- 1 # variance of random effect
#'#-----------------------------------#
#'# Generate function beta(t)-- Linear
#'#-----------------------------------#
#'my.coef <- 10
#'my.beta <- (my.t - min(my.t))/(max(my.t)- min(my.t)) * my.coef
#'#-----------------------------------------------------------------#
#'# Generate read covarage matrix my.X and Z the percentage of type A
#'#-----------------------------------------------------------------#
#'# Generate X (Read Covarage)---using the bootstrap
#'my.X <- matrix(NA, nrow=my.n, ncol = length(pi.0))
#'for ( i in 1:my.n){
#'  my.X[i,] <- sample(y$CT, size = length(pi.0), replace=T)
#'}
#'my.Z <- rbeta(my.n, shape1 = 2, shape2=2)
#'my.data <- BSMethGammSim(n=my.n, posit=my.t, theta.0 = theta.0, beta=my.beta, sigma.ee = sigma.ee,
#'                         X=my.X, Z = my.Z)
#'my.data <- BSMethGammSim(n=my.n, posit=my.t, theta.0 = theta.0, beta=my.beta, sigma.ee = sigma.ee,
#'                         X=my.X, Z = my.Z,random.eff = F)
#'
#' @export
BSMethSim <-
function(n, posit, theta.0, beta, random.eff = F, mu.e=0,
                          sigma.ee=1,p0 = 0.003, p1 = 0.9, X, Z,binom.link="logit"){
  if( !is.matrix((Z)) ) message ("covariate Z is not a matrix")
#  if( !is.matrix(beta) ) message ("the covariate effect parameter beta is not a matrix")

  if( !(nrow(X)==nrow(Z) & nrow(X) == n) ) message("Both X and Z should have n rows")
  if( !(ncol(X)==nrow(theta.0) & ncol(X) ==nrow (beta) & ncol(X)==length(posit) )) message ("The columns of X should be the same as length of beta theta.0 and posit; They all equals to the number of CpGs")
  if( ncol(beta)!= ncol(Z)) message("beta and Z should have the same dimentions")

  # the random effect term
  if(random.eff == T){
    my.e <- rnorm(n, mean=mu.e, sd = sqrt(sigma.ee))
  }else{
    my.e <- rep(mu.e, n)
  }

  my.theta <- t(sapply(1:n, function(i){
    theta.0 + rowSums(sapply(1:ncol(Z), function(j){Z[i,j] * beta[,j]})) + my.e[i]
  }))


  # Transform my.theta to my.pi for each (i, j)

  my.pi <- t(sapply(1:nrow(my.theta), function(i){
    #exp(my.theta[i,])/(1+exp(my.theta[i,]))
    binomial(link=binom.link)$linkinv(my.theta[i,])
  }))
  #  Generate S-ij based on the my.pi and my.Z
  #---------------------------------------------------#
  my.S <- my.pi
  for ( i in 1:nrow(my.S)){
    for ( j in 1:ncol(my.S)){
      my.S[i,j] <- rbinom (1, size = X[i,j], prob= my.pi[i,j])
    }
  }
  #---------------------------------------------------#
  # Generate Y-ij based on the S-ij and the error rate (1-p1) and p0
  #---------------------------------------------------#
  my.Y <- my.S
  for ( i in 1:nrow(my.Y)){
    for ( j in 1:ncol(my.Y)){
      my.Y[i,j] <- sum(rbinom(my.S[i,j], size =1, prob=p1)) +
        sum(rbinom(X[i,j]-my.S[i,j], size = 1, prob=p0))
    }
  }
  out = list(S = my.S, Y = my.Y, theta = my.theta, pi = my.pi)
}
