#' @title A smoothed-EM algorithm to estimate covariate effects and test regional association in Bisulfite Sequencing-derived methylation data
#'
#' @description This function fits a binomial regression model where the outcome - methylated reads- are contaminated by known error rates \code{p0} and \code{p1} and the covariate effects are smoothly varying across genomic positions. The functional parameters of the smooth covariate effects are first represented by a linear combination of a bunch of restricted cubic splines (with dimention \code{n.k}), and a smoothness penalization term which depends on the smoothing parameters \code{lambdas} is also added to control smoothness.
#' @description The estimation is performed by an iterated EM algorithm. Each M step constitutes an outer Newton's iteration to estimate smoothing parameters \code{lambdas} and an inner P-IRLS iteration to estimate spline coefficients \code{alpha} for the covariate effects. Currently, the computation in the M step depends the implementation of \code{gam()} in package \code{mgcv}.
#' @param BEM.obj an output from the function \code{binomRegMethModel}
#' @param newdata the data set whose predictions are calculated
#' @param type return the predicted methylation proportion or the predicted response (in logit or other binom.link scale)
#' @return This function returns the predicted methylation levels
#' @author  Kaiqiong Zhao
#' @examples
#' #------------------------------------------------------------#
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' out <- binomRegMethModel(
#'   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#'   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200, detail=FALSE
#' )
#' plot_binomRegMethModel(out, same.range=FALSE)
#' @export
pred_binomRegMethModel <- function(BEM.obj, newdata=NULL, type="proportion") {
    uni.pos <- BEM.obj$uni.pos
    covs <- colnames(BEM.obj$Beta.out)
    if (!type %in% c("proportion", "link.scale")) {
        stop("type should be either proportion or link.scale")
    }

    if (is.null(newdata)) {
        if (type == "proportion") {
            return(BEM.obj$est.pi)
        }
        if (type == "link.scale") {
            return(log(BEM.obj$est.pi/(1 - BEM.obj$est.pi)))
        }
    } else {
        newdata <- data.frame(newdata, Intercept=1)
        if (all(setdiff(colnames(newdata), "Position") %in% colnames(BEM.obj$Beta.out))) {
            id <- match(newdata$Position, uni.pos)

            if (any(is.na(id))) {
                stop("The positions in the newdata should be the exactly the same as the positions fited in object BEM.obj")
            }
            ## estimated beta(t) for each t in the role of newdata
            beta.s <- BEM.obj$Beta.out[id, ]

            newdata <- newdata[, -which(colnames(newdata) == "Position")]
            newdata <- newdata[, match(covs, colnames(newdata))]

            pred.link <- sapply(seq_len(nrow(newdata)), function(i) {
                sum(beta.s[i, ] * newdata[i, ])
            })

            if (type == "link.scale") {
                return(pred.link)
            }
            if (type == "proportion") {
                pred.pi <- exp(pred.link)/(1 + exp(pred.link))
                return(pred.pi)
            }
        } else {
            stop("The covariates used to fit object BEM.obj should appeared in newdata")
        }
    }
}
