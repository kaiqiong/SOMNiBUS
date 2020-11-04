#' @title Plot the smooth covariate effect
#'
#' @description This function accepts an output object from function
#' \code{binomRegMethModel} and print out a plot of the estimated
#' covariate effect
#' across the region for each test covariate.
#' @param BEM.obj an output object from function
#' \code{binomRegMethModel}
#' @param mfrow the plot parameters to specify the layout of each plot
#' @param same.range specify whether the plots should be in the same
#' vertical scale
#' @return This function prints out a plot of smooth covariate effects
#' and
#' its pointwise confidence intervals
#' @author  Kaiqiong Zhao
#' @importFrom graphics abline axis lines par plot
#' @examples
#' #------------------------------------------------------------#
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' out <- binomRegMethModel(
#'   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#'   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200, detail=FALSE
#' )
#' binomRegMethModelPlot(out, same.range=FALSE)
#' @export
binomRegMethModelPlot <- function(BEM.obj, mfrow = NULL, same.range = FALSE) {
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
        yylim <- t(vapply(seq_len(ncovs), function(i) {
            yylim[i, ] <- c(ifelse(min(ll[, i]) > 0, 0, min(ll[, i])), 
                max(hh[, i]))
        }, rep(1, 2)))
    }
    
    for (ii in seq_len(ncovs)) {
        plot(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], BEM.obj$Beta.out[order(BEM.obj$uni.pos), 
            ii], col = "red", xaxt = "n", type = "l", xlab = "Genomic Position", 
            ylab = " ", main = covs.names[ii], lwd = 2, ylim = yylim[ii, 
                ])
        axis(side = 1, at = BEM.obj$uni.pos[order(BEM.obj$uni.pos)], labels = FALSE, 
            lwd = 0.5, lwd.ticks = 0.5, tck = 0.03)
        axis(side = 1, at = seq(round(min(BEM.obj$uni.pos)), round(max(BEM.obj$uni.pos)), 
            length.out = 10), tck = -0.02)
        lines(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], ll[order(BEM.obj$uni.pos), 
            ii], lty = 2, col = "red")
        lines(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], hh[order(BEM.obj$uni.pos), 
            ii], lty = 2, col = "red")
        abline(h = 0, lty = 2)
    }
}
