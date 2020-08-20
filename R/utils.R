#' @title hessianComp compute the Hessian matrix
#'
#' @description compute the Hessian matrix for an EM estimate
#' @param w_ij the diagnal values of the weight matrix
#' @param new.par estimate of alpha
#' @param new.lambda estimate of lambda
#' @param X a vector of read depths
#' @param Y a vector of methylated counts
#' @param my.design.matrix design matrix from the final fit
#' @param gam_smoothMat the smooth matrix from the final gam fit
#' @param Z covariate matrix
#' @param pred.pi predicted methylation probability from the final fit
#' @param p0 the probability of observing a methylated read when the
#' underlying
#' true status is unmethylated. \code{p0} is the rate
#' of false methylation calls, i.e. false positive rate.
#' @param p1 the probability of observing a methylated read when the
#' underlying
#' true status is methylated. \code{1-p1} is the rate
#' of false non-methylation calls, i.e. false negative rate.
#' @param disp_est estimated dispersion parameter
#' @param RanEff whether a subject-level Random effect is added or not
#' @param N number of unique samples in the provided dataset
#' @return a Hessian matrix for all alphas (including the alphas for RE)
#' @references Zhao, Kaiqiong, et al. 'A novel statistical method for
#' modeling  covariate effects in bisulfite sequencing derived measures
#' of DNA methylation.'Biometrics (2020).
#' @author  Kaiqiong Zhao
#' @importFrom Matrix bdiag
#' @noRd
hessianComp <- function(w_ij, new.par, new.lambda, X, Y, my.design.matrix, 
    gam_smoothMat, Z, pred.pi, p0, p1, disp_est, RanEff, N) {
    
    ## Q1: the second partial derivative w.r.t alpha^2 Q2: the second
    ## derivative w.r.t alpha & alpha_star
    res <- outer(seq_len(length(new.par)), seq_len(length(new.par)), Vectorize(function(l, 
        m) {
        sum(-X * w_ij * my.design.matrix[, m] * my.design.matrix[, l])
    }))
    smoth.mat <- lapply(as.list(seq_len(ncol(Z) + 1)), function(i) {
        gam_smoothMat[[i]]$S[[1]] * new.lambda[i]
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
