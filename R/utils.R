#' @title hessianComp Lorem ipsum dolor sit amet
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
#' @author  Kaiqiong Zhao
#' @importFrom Matrix bdiag
#' @noRd
hessianComp <- function(w_ij, new.par, new.lambda, X, Y, my.design.matrix,
    gam.int, Z, pred.pi, p0, p1, disp_est, RanEff, N) {

    ## Q1: the second partial derivative w.r.t alpha^2 Q2: the second
    ## derivative w.r.t alpha & alpha_star
    res <- outer(seq_len(length(new.par)), seq_len(length(new.par)), Vectorize(function(l,
        m) {
        sum(-X * w_ij * my.design.matrix[, m] * my.design.matrix[, l])
    }))
    smoth.mat <- lapply(as.list(seq_len(ncol(Z) + 1)), function(i) {
        gam.int$smooth[[i]]$S[[1]] * new.lambda[i]
    })  ## extract the penalty matrix
    ## assume the lambda for the constant of the intercept is 0 -- no
    ## penalization
    smoth.mat[[length(smoth.mat) + 1]] <- 0
    if (RanEff) {
        ## !!!! Otherwise, we get very wide CI
        smoth.mat[[length(smoth.mat) + 1]] <- diag(N) * new.lambda[ncol(Z) +
            2]
        span.penal.matrix <- as.matrix(Matrix::bdiag(smoth.mat[c(length(smoth.mat) -
            1, (seq_len((length(smoth.mat) - 2))), length(smoth.mat))]))
    } else {
        span.penal.matrix <- as.matrix(Matrix::bdiag(smoth.mat[c(length(smoth.mat),
            (seq_len((length(smoth.mat) - 1))))]))
    }

    Q1_with_lambda <- res - span.penal.matrix/disp_est
    Q1_no_lambda <- res

    Q2 <- outer(seq_len(length(new.par)), seq_len(length(new.par)), Vectorize(function(l,
        m) {
        term1 <- Y * p1 * p0/(p1 * pred.pi + p0 * (1 - pred.pi))^2 + (X -
            Y) * (1 - p1) * (1 - p0)/((1 - p1) * pred.pi + (1 - p0) * (1 -
            pred.pi))^2
        sum(term1 * w_ij * my.design.matrix[, m] * my.design.matrix[, l])
    }))

    return(Q1_with_lambda + Q2)
}
