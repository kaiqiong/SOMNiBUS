#' @title A smoothed-EM algorithm to estimate covariate effects and test regional association in Bisulfite Sequencing-derived methylation data
#'
#' @description This function fits a binomial regression model where the outcome - methylated reads- are contaminated by known error rates \code{p0} and \code{p1} and the covariate effects are smoothly varying across genomic positions. The functional parameters of the smooth covariate effects are first represented by a linear combination of a bunch of restricted cubic splines (with dimention \code{n.k}), and a smoothness penalization term which depends on the smoothing parameters \code{lambdas} is also added to control smoothness.
#' @description The estimation is performed by an iterated EM algorithm. Each M step constitutes an outer Newton's iteration to estimate smoothing parameters \code{lambdas} and an inner P-IRLS iteration to estimate spline coefficients \code{alpha} for the covariate effects. Currently, the computation in the M step depends the implementation of \code{gam()} in package \code{mgcv}.
#' @param BEM.obj an output object from function \code{BSMethEM}
#' @param mfrow the plot parameters to specify the layout of each plot
#' @param same.range specify whether the plots should be in the same vertical scale
#' @return This function return a plot of smooth covariate effects and its pointwise confidence intervals
#' @author  Kaiqiong Zhao
#' @importFrom graphics abline axis lines par plot
#' @examples
#' #------------------------------------------------------------#
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' out <- BSMethEM(
#'   data = RAdat.f, n.k = rep(5, 3), p0 = 0.003307034, p1 = 0.9,
#'   epsilon = 10^(-6), epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE
#' )
#' plot_BSMethEM(out, same.range = FALSE)
#' @export
plot_BSMethEM <- function(BEM.obj, mfrow = NULL, same.range = FALSE) {
  ncovs <- ncol(BEM.obj$Beta.out)

  if (is.null(mfrow)) {
    par(mfrow = c(1, ncovs), mar = c(4, 2.6, 3, 1))
  } else {
    par(mfrow = mfrow, mar = c(4, 2.6, 3, 1))
  }
  covs.names <- rownames(BEM.obj$reg.out)
  ll <- BEM.obj$Beta.out - 1.96 * BEM.obj$SE.out
  hh <- BEM.obj$Beta.out + 1.96 * BEM.obj$SE.out

  yylim <- matrix(NA, nrow = ncovs, ncol = 2)
  if (same.range) {
    yylim <- matrix(rep(c(min(ll), max(hh)), ncovs), ncol = 2, byrow = TRUE)
  } else {
    yylim <- t(sapply(1:ncovs, function(i) {
      yylim[i, ] <- c(ifelse(min(ll[, i]) > 0, 0, min(ll[, i])), max(hh[, i]))
    }))
  }

  for (ii in 1:ncovs) {
    plot(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], BEM.obj$Beta.out[
      order(BEM.obj$uni.pos),
      ii
    ],
    col = "red", xaxt = "n", type = "l", xlab = "Genomic Position",
    ylab = " ", main = covs.names[ii], lwd = 2, ylim = yylim[ii, ]
    )
    axis(
      side = 1, at = BEM.obj$uni.pos[order(BEM.obj$uni.pos)], labels = FALSE,
      lwd = 0.5, lwd.ticks = 0.5, tck = 0.03
    )
    axis(side = 1, at = seq(round(min(BEM.obj$uni.pos)), round(max(BEM.obj$uni.pos)),
      length.out = 10
    ), tck = -0.02)
    lines(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], ll[
      order(BEM.obj$uni.pos),
      ii
    ], lty = 2, col = "red")
    lines(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], hh[
      order(BEM.obj$uni.pos),
      ii
    ], lty = 2, col = "red")
    abline(h = 0, lty = 2)
  }
}
