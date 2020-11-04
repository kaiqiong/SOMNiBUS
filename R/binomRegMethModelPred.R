#' @title A smoothed-EM algorithm to estimate covariate effects and test
#' regional association in Bisulfite Sequencing-derived methylation data
#'
#' @description This function returns the predicted methylation levels
#' @param BEM.obj an output from the function \code{binomRegMethModel}
#' @param newdata the data set whose predictions are calculated; with
#' columns
#' 'Position', and covariate names that can be matched to the BEM.obj
#' @param type return the predicted methylation proportion or the
#' predicted
#' response (in logit or other binom.link scale)
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
#' binomRegMethModelPred(out)
#' @export
binomRegMethModelPred <- function(BEM.obj, newdata = NULL, type = "proportion") {
    uni.pos <- BEM.obj$uni.pos
    covs <- colnames(BEM.obj$Beta.out)
    
    binomRegMethModelPredChecks(type = type)
    
    if (is.null(newdata)) {
        if (type == "proportion") {
            return(BEM.obj$est.pi)
        }
        if (type == "link.scale") {
            return(log(BEM.obj$est.pi/(1 - BEM.obj$est.pi)))
        }
    } else {
        newdata <- data.frame(newdata, Intercept = 1)
        if (all(colnames(BEM.obj$Beta.out) %in% colnames(newdata))) {
            id <- match(newdata$Position, uni.pos)
            
            if (any(is.na(id))) {
                stop("The positions in the newdata should be the exactly
                     the same as the positions fited in object BEM.obj")
            }
            ## estimated beta(t) for each t in the role of newdata
            beta.s <- BEM.obj$Beta.out[id, ]
            
            newdata <- newdata[, -which(colnames(newdata) == "Position")]
            newdata <- newdata[, match(covs, colnames(newdata))]
            
            pred.link <- vapply(seq_len(nrow(newdata)), function(i) {
                sum(beta.s[i, ] * newdata[i, ])
            }, 1)
            
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

#' @title Some checks for binomRegMethModelPred
#'
#' @description check whether an appropriate prediction type
#' is given or not
#' @param type return the predicted methylation proportion or
#' the predicted response (in logit or other binom.link scale)
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
binomRegMethModelPredChecks <- function(type) {
    if (!type %in% c("proportion", "link.scale")) {
        stop("type should be either proportion or link.scale")
    }
}
