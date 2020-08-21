#' @title A smoothed-EM algorithm to estimate covariate effects and test
#' regional association in Bisulfite Sequencing-derived methylation data
#' @description This function fits a (dispersion-adjusted) binomial
#' regression model to regional methylation data, and reports the
#' estimated smooth covariate effects and regional p-values for the
#' test of DMRs (differentially methylation regions). Over or under
#' dispersion across loci is accounted for in the model by the combination
#' of a multiplicative dispersion parameter (or scale parameter) and a
#' sample-specific random effect.
#' @description This method can deal with outcomes, i.e. the number of
#' methylated reads in a region, that are contaminated by known
#' false methylation calling rate (\code{p0}) and false non-methylation
#' calling rate (\code{1-p1}).
#' @description The covariate effects are assumed to smoothly vary across
#' genomic regions. In order to estimate them, the algorithm first
#' represents the functionalparamters by a linear combination
#' of a set of restricted cubic splines (with dimention
#' \code{n.k}), and a smoothness penalization term
#' which depends on the smoothing parameters \code{lambdas} is also
#' added to control smoothness.
#' The estimation is performed by an iterated EM algorithm. Each
#' M step constitutes
#' an outer Newton's iteration to estimate smoothing parameters
#' \code{lambdas} and an
#' inner P-IRLS iteration to estimate spline coefficients
#' \code{alpha} for the covariate
#' effects. Currently, the computation in the M step depends the
#' implementation of
#' \code{gam()} in package \code{mgcv}.
#' @param data a data frame with rows as individual CpGs appeared
#' in all the samples. The
#' first 4 columns should contain the information of `Meth_Counts`
#' (methylated counts),
#' `Total_Counts` (read depths), `Position` (Genomic position for
#' the CpG site) and `ID`
#' (sample ID). The covariate information, such as disease status
#' or cell type composition,
#' are listed in column 5 and onwards.
#' @param n.k a vector of basis dimensions for the intercept and
#' individual covariates.
#' \code{n.k} specifies an upper limit of the degrees of each
#' functional parameters.
#' @param p0 the probability of observing a methylated read when
#' the underlying true
#' status is unmethylated. \code{p0} is the rate of false
#' methylation calls, i.e.
#' false positive rate.
#' @param p1 the probability of observing a methylated read when
#' the underlying true
#' status is methylated. \code{1-p1} is the rate of false
#' non-methylation calls, i.e.
#' false negative rate.
#' @param covs a vector of covariate names. The covariates with
#' names in \code{covs} will
#' be included in the model and their covariate effects will be
#' estimated. The default
#' is to fit all covariates in \code{data}
#' @param Quasi whether a Quasi-likelihood estimation approach
#' will be used; in other
#' words, whether a multiplicative dispersion is added in the
#' model or not.
#' @param epsilon numeric; stopping criterion for the closeness of
#' estimates of spline
#' coefficients from two consecutive iterations.
#' @param epsilon.lambda numeric; stopping criterion for the
#' closeness of estimates of smoothing parameter \code{lambda}
#' from two consecutive iterations.
#' @param maxStep the algorithm will step if the iteration steps
#' exceed \code{maxStep}
#' @param detail indicate whether print the number of iterations
#' @param binom.link the link function used in the binomial
#' regression model; the default
#' is the logit link
#' @param method the method used to estimate the smoothing
#' parameters. The default is the
#' 'REML' method which is generally better than prediction based
#' criterion \code{GCV.cp}
#' @param RanEff whether sample-level random effects are added or not
#' @param reml.scale whether a REML-based scale (dispersion)
#' estimator is used. The
#' default is Fletcher-based estimator
#' @param scale nagative values mean scale paramter should be
#' estimated; if a positive
#' value is provided, a fixed scale will be used.
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{est}: estimates of the spline basis coefficients
#' \code{alpha}
#' \item \code{lambda}: estimates of the smoothing parameters
#' for each functional paramters
#' \item \code{est.pi}: predicted methylation levels for each
#' row in the input \code{data}
#' \item \code{ite.points}: estimates of \code{est}, \code{lambda}
#' at each EM iteration
#' \item \code{cov1}: estimated variance-covariance matrix of the
#' basis coefficients \code{alphas}
#' \item \code{reg.out}: regional testing output obtained using
#' Fletcher-based dispersion
#' estimate; an additional 'ID' row would appear if RanEff is TRUE
#' \item \code{reg.out.reml.scale}:regional testing output obtained
#' sing REML-based
#' dispersion estimate;
#' \item \code{reg.out.gam}:regional testing output obtained using
#' (Fletcher-based)
#' dispersion estimate from mgcv package;
#' \item \code{phi_fletcher}: Fletcher-based estimate of the
#' (multiplicative) dispersion
#' parameter
#' \item \code{phi_reml}: REML-based estimate of the (multiplicative)
#' dispersion parameter
#' \item \code{phi_gam}: Estimated dispersion parameter reported by
#' mgcv
#' \item \code{SE.out}: a matrix of the estimated pointwise Standard
#' Errors (SE); number
#' of rows are the number of unique CpG sites in the input data and
#' the number of columns
#' equal to the total number of covariates fitted in the model
#' (the first one is the intercept)
#' \item \code{SE.out.REML.scale}: a matrix of the estimated
#' pointwise Standard Errors (SE);
#' the SE calculated from the REML-based dispersion estimates
#' \item \code{uni.pos}: the genomic postions for each row of
#' CpG sites in the matrix \code{SE.out}
#' \item \code{Beta.out}: a matrix of the estimated covariate
#' effects beta(t), here t
#' denots the genomic positions.
#' \item \code{ncovs}: number of functional paramters in the model
#' (including the
#' intercept)
#' \item \code{sigma00}: estimated variance for the random effect
#' if RanEff is TRUE;
#' NA if RanEff is FALSE
#' }
#' @author  Kaiqiong Zhao
#' @seealso  \link[mgcv]{gam}
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' out <- binomRegMethModel(
#'   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#'   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#'   detail=FALSE
#' )
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom mgcv model.matrix.gam
#' @importFrom mgcv s
#' @importFrom Matrix bdiag
#' @importFrom stats as.formula binomial pchisq rbinom rnorm
#' quasibinomial residuals predict model.matrix runif
#'
#' @export
binomRegMethModel <- function(data, n.k, p0 = 0.003, p1 = 0.9, Quasi = TRUE,
    epsilon = 10^(-6), epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE,
    binom.link = "logit", method = "REML", covs = NULL, RanEff = TRUE,
    reml.scale = FALSE, scale = -2) {
    initOut <- binomRegMethModelInit(data = data, covs = covs, n.k = n.k)
    Z <- initOut$Z
    fitGamOut <- fitGam(data = initOut$data, Quasi = Quasi, binom.link = binom.link,
        method = method, RanEff = RanEff, scale = scale, Z = Z, n.k = n.k)
    phi_fletcher <- phiFletcher(data = fitGamOut$data, Quasi = Quasi, reml.scale = reml.scale,
        scale = scale, gam.int = fitGamOut$gam.int)
    fitEMOut <- fitEM(p0=p0, p1=p1, fitGamOut=fitGamOut, n.k=n.k, binom.link=binom.link, method=method, Z=Z, Quasi=Quasi,
        scale=scale, reml.scale=reml.scale, phi_fletcher=phi_fletcher, epsilon=epsilon, maxStep=maxStep, detail=detail)
    out <- fitEMOut$out
    Est.points <- fitEMOut$Est.points
    phi_fletcher <- out$phi_fletcher
    my.design.matrix <- mgcv::model.matrix.gam(out$GamObj)
    lengthUniqueDataID <- length(unique(fitGamOut$data$ID))
    phi_reml <- out$GamObj$reml.scale
    estimateBZOut <- estimateBZ(Posit = fitGamOut$data$Posit, my.design.matrix = my.design.matrix,
        ncolsZ = ncol(Z), n.k = n.k)
    Beta.out <- estimateBeta(BZ = estimateBZOut$BZ, BZ.beta = estimateBZOut$BZ.beta,
        n.k = n.k, Z = Z, out = out)
    estimateVarOut <- estimateVar(out=out, phi_fletcher=phi_fletcher, fitGamOut=fitGamOut, my.design.matrix=my.design.matrix,
        Z=Z, p0=p0, p1=p1, RanEff=RanEff, lengthUniqueDataID=lengthUniqueDataID, n.k=n.k)
    estimateSEOut <- estimateSE(estimateBZOut=estimateBZOut, Z=Z, estimateVarOut=estimateVarOut, phi_fletcher=phi_fletcher,
        phi_reml=phi_reml)
    X_d <- extractDesignMatrix(GamObj = out$GamObj)
    edf1.out <- out$edf1  ## and p value calculation tr(2A - A^2)
    edf.out <- out$edf  ## Effective degree of freedom: edf --trace of the hat matrix
    resi_df <- nrow(fitGamOut$data) - sum(edf.out)
    s.table <- binomRegMethModelSummary(GamObj=out$GamObj, var.cov.alpha=estimateVarOut$var.cov.alpha,
        new.par=out$par, edf.out=edf.out, edf1.out=edf1.out, X_d=X_d, resi_df=resi_df, Quasi=Quasi, scale=scale, RanEff=RanEff,
        Z=Z)
    s.table.REML.scale <- binomRegMethModelSummary(GamObj=out$GamObj, var.cov.alpha=estimateVarOut$var.cov.alpha/phi_fletcher *
        phi_reml, new.par=out$par, edf.out=edf.out, edf1.out=edf1.out, X_d=X_d, resi_df=resi_df, Quasi=Quasi, scale=scale,
        RanEff=RanEff, Z=Z)
    if (RanEff) {
        sigma00 <- out$GamObj$reml.scale/out$GamObj$sp["s(ID)"]
    } else {
        sigma00 <- NA
    }
    return(out <- list(est = out$par, lambda = out$lambda, est.pi = out$pi.ij,
        Beta.out = Beta.out, phi_fletcher = phi_fletcher, phi_reml = phi_reml,
        phi_gam = out$GamObj$scale, reg.out = s.table, reg.out.reml.scale = s.table.REML.scale,
        cov1 = estimateVarOut$var.cov.alpha, reg.out.gam = summary(out$GamObj)$s.table,
        SE.out = estimateSEOut$SE.out, SE.out.REML.scale = estimateSEOut$SE.out.REML.scale,
        uni.pos = estimateBZOut$uni.pos, ncovs = ncol(Z) + 1, ite.points = Est.points,
        sigma00 = sigma00))
}
#' @title Run EM algorithm to obtain the estimate of alpha
#'
#' @description Run EM algorithm to obtain the estimate of alpha
#' @param p0 false positive error rates
#' @param p1 1-p1: false negative error rates
#' @param fitGamOut an output from fitGam
#' @param n.k number of knots for all covariates (including intercept)
#' @param binom.link link for binomial GLM
#' @param method method to estimate lambda
#' @param Z covariate matrix
#' @param Quasi whether quasibinomial is used
#' @param scale negative means unknown scale; positive means fixed scale
#' @param reml.scale whether a REML-based scale estimate is used
#' @param phi_fletcher Fletcher-based scale estimate
#' @param epsilon stopping criteria
#' @param maxStep maximum iteration steps
#' @param detail whether to print out the iteration
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{out} a list
#' \item \code{Est.points} a matrix storing the estimate in each step
#' }
#' @author Kaiqiong Zhao
#' @noRd
fitEM <- function(p0, p1, fitGamOut, n.k, binom.link, method, Z, Quasi,
    scale, reml.scale, phi_fletcher, epsilon, maxStep, detail) {
    if (p0 > 0 | p1 < 1) {
        out <- binomRegMethModelUpdate(data = fitGamOut$data, pi.ij = fitGamOut$gam.int$fitted.values,
            p0 = p0, p1 = p1, n.k = n.k, binom.link = binom.link, method = method,
            Z = Z, my.covar.fm = fitGamOut$my.covar.fm, Quasi = Quasi,
            scale = scale, reml.scale = reml.scale)
        Est.points <- rbind(c(fitGamOut$gam.int$coefficients, fitGamOut$gam.int$sp,
            phi_fletcher), c(out$par, out$lambda, out$phi_fletcher))
        old.pi.ij <- fitGamOut$gam.int$fitted.values
        iter <- 1
        while (sqrt(sum((out$pi.ij - old.pi.ij)^2)) > epsilon & iter <
            maxStep) {
            if (detail) {
                message(paste0("iteration", iter))
            }
            iter <- iter + 1
            old.pi.ij <- out$pi.ij
            out <- binomRegMethModelUpdate(data = fitGamOut$data, pi.ij = old.pi.ij,
                p0 = p0, p1 = p1, n.k = n.k, binom.link = binom.link, method = method,
                Z = Z, my.covar.fm = fitGamOut$my.covar.fm, Quasi = Quasi,
                scale = phi_fletcher)
            Est.points <- rbind(Est.points, c(out$par, out$lambda, out$phi_fletcher))
        }
    } else {
        out <- list(pi.ij = fitGamOut$gam.int$fitted.values, par = fitGamOut$gam.int$coefficients,
            lambda = fitGamOut$gam.int$sp, edf1 = fitGamOut$gam.int$edf1,
            pearson_res = residuals(fitGamOut$gam.int, type = "pearson"),
            deviance_res = residuals(fitGamOut$gam.int, type = "deviance"),
            edf = fitGamOut$gam.int$edf, phi_fletcher = phi_fletcher, GamObj = fitGamOut$gam.int)
        Est.points <- c(fitGamOut$gam.int$coefficients, fitGamOut$gam.int$sp,
            phi_fletcher)
    }
    return(list(out = out, Est.points = Est.points))
}

#' @title Get the basis matrix for beta0(t) and beta(t)
#'
#' @description This function aims to extract the basis matrix BZ,
#' and BZ.beta from the expansion beta0(t) = BZ * alpha0
#' and beta(t) = BZ.beta * alpha.
#' @param Posit  the column 'Posit' from the input data frame
#' \code{data}
#' @param my.design.matrix the design matrix for a \code{fitGamOut}
#' with data as input
#' and Z as covariates.
#' @param ncolsZ number of columns of covariate matrix Z
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{uni.pos} the unique genomic positions present in
#' the data
#' \item \code{BZ} basis matrix for the intercept beta0(t); a
#' matrix with
#' \code{length(uni.pos)} rows and \code{length(alpha0)} columns
#' \item \code{BZ.beta} a list of length \code{ncolsZ}; each
#' corresponds to the basis
#' matrix for beta_p(t), p = 1, 2, .. \code{ncolsZ}
#' }
#' @author  Kaiqiong Zhao, Simon Laurin-Lemay
#' @import mgcv
#' @noRd
estimateBZ <- function(Posit, my.design.matrix, ncolsZ, n.k) {
    uni.pos <- unique(Posit)
    uni.id <- match(uni.pos, Posit)
    BZ <- my.design.matrix[uni.id, seq_len(n.k[1])]
    BZ.beta <- lapply(seq_len(ncolsZ), function(i) {
        mgcv::smooth.construct(mgcv::s(Posit, k = n.k[i + 1], fx = FALSE,
            bs = "cr"), data = data.frame(Posit = Posit[uni.id]), knots = NULL)$X
    })
    return(out <- list(uni.pos = uni.pos, BZ = BZ, BZ.beta = BZ.beta))
}

#' @title estimate of beta(t)
#'
#' @description Given a final GAM output, extract the estimates of
#' beta(t) = BZ * alpha, where t is the unique genomic positions
#' present in the input data; only extract the fixed effect estimates
#' @param BZ basis matrix for the intercept beta0(t); a matrix with
#' \code{length(uni.pos)} rows and \code{length(alpha0)} columns
#' @param BZ.beta a list of length \code{ncolsZ}; each corresponds
#' to the basis matrix
#' for beta_p(t)
#' @param n.k a vector of basis dimensions for the intercept and
#' individual
#' covariates. \code{n.k} specifies an upper limit of the degrees
#' of each functional
#' parameters.
#' @param Z a covariate matrix
#' @param out an object from \code{binomRegMethModelUpdate}
#' @return \code{Beta.out}: a matrix of the estimated covariate
#' effects beta(t),
#' here t denots the unique genomic positions for intercept
#' and all covariates. Beta.out = cbind(beta.0(t), beta.1(t), ...
#' beta.P(t));
#' \code{length(uni.pos)} rows and \code{ncol(Z)+1} columns
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
estimateBeta <- function(BZ, BZ.beta, n.k, Z, out) {
    ## ------ 3: estimate of beta(t) --- #
    cum_s <- cumsum(n.k)
    alpha.sep <- lapply(seq_len(ncol(Z)), function(i) {
        out$par[(cum_s[i] + 1):cum_s[i + 1]]
    })
    alpha.0 <- out$par[seq_len(n.k[1])]

    Beta.out <- cbind(BZ %*% alpha.0, vapply(seq_len(ncol(Z)), function(i) {
        BZ.beta[[i]] %*% alpha.sep[[i]]
    }, FUN.VALUE = rep(1, nrow(BZ))))

    colnames(Beta.out) <- c("Intercept", colnames(Z))
    return(Beta.out)
}
#' @title Estimate the variance covariance matrix
#'
#' @description Estimate the variance covariance matrix
#' @param out an output from the update function
#' @param phi_fletcher Fletcher-based scale estimate
#' @param fitGamOut an output from fitGam
#' @param RanEff whether a RE is considered
#' @param n.k number of knots for all covariates (including intercept)
#' @param lengthUniqueDataID number of samples in the data
#' @param p0 false positive error rates
#' @param p1 1-p1: false negative error rates
#' @param Z covariate matrix
#' @param my.design.matrix design matrix for the all data
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{var.cov.alpha} var of alpha
#' \item \code{var.alpha.0} var of alpha0
#' \item \code{var.alpha.sep} var of alpha_p, p = 1, 2, P
#' }
#' @author Kaiqiong Zhao
#' @noRd
estimateVar <- function(out, phi_fletcher, fitGamOut, my.design.matrix,
    Z, p0, p1, RanEff, lengthUniqueDataID, n.k) {
    cum_s <- cumsum(n.k)
    w_ij <- out$pi.ij * (1 - out$pi.ij)/phi_fletcher
    H <- hessianComp(w_ij = w_ij, new.par = out$par, new.lambda = out$lambda,
        X = fitGamOut$data$X, Y = fitGamOut$data$Y, my.design.matrix = my.design.matrix,
        gam_smoothMat = out$GamObj$smooth, Z = Z, pred.pi = out$pi.ij,
        p0 = p0, p1 = p1, disp_est = phi_fletcher, RanEff = RanEff, N = lengthUniqueDataID)

    var.cov.alpha <- solve(-H)
    var.alpha.0 <- var.cov.alpha[seq_len(n.k[1]), seq_len(n.k[1])]
    var.alpha.sep <- lapply(seq_len(ncol(Z)), function(i) {
        var.cov.alpha[(cum_s[i] + 1):cum_s[i + 1], (cum_s[i] + 1):cum_s[i +
            1]]
    })
    return(list(var.alpha.0 = var.alpha.0, var.alpha.sep = var.alpha.sep,
        var.cov.alpha = var.cov.alpha))
}
#' @title estimate SE of beta(t) for each covariates
#'
#' @description estimate SE of beta(t) for each covariates
#' @param estimateBZOut a final Gam Object
#' @param Z the ith smooth/covariate to be evaluated
#' @param estimateVarOut output from estimateVar
#' @param phi_fletcher Fletcher-based dispersion estimate
#' @param phi_reml REML-based dispersion estimate
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{SE.out} SE based on phi_fletcher
#' \item \code{SE.out.REML.scale} SE based on phi_reml
#' }
#' @author Kaiqiong Zhao
#' @noRd
estimateSE <- function(estimateBZOut, Z, estimateVarOut, phi_fletcher,
    phi_reml) {
    var.alpha.0 <- estimateVarOut$var.alpha.0
    var.alpha.sep <- estimateVarOut$var.alpha.sep
    SE.out <- cbind(sqrt(pmax(0, rowSums((estimateBZOut$BZ %*% var.alpha.0) *
        estimateBZOut$BZ))), vapply(seq_len(ncol(Z)), function(i) {
        sqrt(pmax(0, rowSums((estimateBZOut$BZ.beta[[i]] %*% var.alpha.sep[[i]]) *
            estimateBZOut$BZ.beta[[i]])))
    }, FUN.VALUE = rep(1, length(estimateBZOut$uni.pos))))

    rownames(SE.out) <- estimateBZOut$uni.pos
    colnames(SE.out) <- c("Intercept", colnames(Z))
    SE.out.REML.scale <- SE.out/sqrt(phi_fletcher) * sqrt(phi_reml)
    return(list(SE.out = SE.out, SE.out.REML.scale = SE.out.REML.scale))
}
#' @title binomRegMethModelSummary Extract the regional testing
#' results from a GamObj
#'
#' @description Extract the regional testing results from a
#' GamObj
#' @param GamObj an fit from mgcv::gam
#' @param var.cov.alpha Estimated variance-covariance matrix
#' for coefficients alpha
#' @param new.par Estimated coefficients alpha
#' @param edf.out Effective degrees of freedom for each alpha,
#' calculated from the
#' trace of the hat matrix
#' @param edf1.out Effective degrees of freedom 2, calculated
#' from tr(2A - A^2); good
#' for chisquare test
#' @param X_d R from the QR decomposition for the design matrix
#' @param resi_df Residual degrees of freedom; nrow(data) -
#' sum(edf.out)
#' @param Quasi Whether a multiplicative dispersion parameter
#' is added or not;
#' quasibinomial or binomial
#' @param scale nagative values mean scale paramter should be
#' estimated;
#' if a positive value is provided, a fixed scale will be used.
#' @param RanEff Whether a subject random effect is added or not
#' @param Z covariate matrix; \code{Z = data[,-c(1:4)]}
#' @return s.table a table for each covariate
#' @seealso  \link[mgcv]{summary.gam}
#' @import mgcv
#' @noRd
binomRegMethModelSummary <- function(GamObj, var.cov.alpha, new.par, edf.out,
    edf1.out, X_d, resi_df, Quasi, scale, RanEff, Z) {
    ii <- 0
    m <- length(GamObj$smooth)
    df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
    for (i in seq_len(m)) {
        outOneSmooth <- regResOneSmooth(GamObj=GamObj, i=i, RanEff=RanEff, Quasi=Quasi, resi_df=resi_df,
            var.cov.alpha=var.cov.alpha, new.par=new.par, edf.out=edf.out, edf1.out=edf1.out, X_d=X_d)
        res <- outOneSmooth$res
        edf1i <- outOneSmooth$edf1i
        edfi <- outOneSmooth$edfi
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
            df <- df[seq_len(ii)]
            chi.sq <- chi.sq[seq_len(ii)]
            edf1 <- edf1[seq_len(ii)]
            edf <- edf[seq_len(ii)]
            s.pv <- s.pv[seq_len(ii)]
        }
        if (!Quasi || scale >= 0) {
            s.table <- cbind(edf, df, chi.sq, s.pv)
            dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df",
                "Chi.sq", "p-value"))
        } else {
            s.table <- cbind(edf, df, chi.sq/df, s.pv)
            dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df",
                "F", "p-value"))
        }
    }
    s.table <- s.table[, -1]
    if (RanEff) {
        rownames(s.table) <- c("Intercept", colnames(Z), "ID")
    } else {
        rownames(s.table) <- c("Intercept", colnames(Z))
    }
    colnames(s.table)[1] <- "EDF"
    return(s.table)
}

#' @title Obtain res for one smooth
#'
#' @description Get the regional summary results for an individual
#' covariate effect
#' @param GamObj a final Gam Object
#' @param i the ith smooth/covariate to be evaluated
#' @param RanEff whether a subject-level RE is added or not
#' @param Quasi whether a quasibinomial dist. is used
#' @param resi_df residual degrees of freedom
#' @param X_d the design matrix (R) for all the smooth
#' @param var.cov.alpha the variance covariance matrix for all the smooth
#' @param new.par the coefficients all the smooth
#' @param edf.out edf1 for all the smooth
#' @param edf1.out edf1 for all the smooth
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{res} output from mgcv:::testStat or
#' mgcv:::reTest
#' \item \code{edfi} edf for the ith smooth
#' \item \code{edf1i} edf1 for the ith smooth
#' }
#' @author Kaiqiong Zhao
#' @noRd
regResOneSmooth <- function(GamObj, i, RanEff, Quasi, resi_df, var.cov.alpha,
    new.par, edf.out, edf1.out, X_d) {
    start <- GamObj$smooth[[i]]$first.para
    stop <- GamObj$smooth[[i]]$last.para
    V <- var.cov.alpha[start:stop, start:stop, drop = FALSE]  ## Bayesian
    p <- new.par[start:stop]  ## params for smooth
    edfi <- sum(edf.out[start:stop])  ## edf for this smooth
    edf1i <- sum(edf1.out[start:stop])
    Xt <- X_d[, start:stop, drop = FALSE]
    fx <- if (inherits(GamObj$smooth[[i]], "tensor.smooth") && !is.null(GamObj$smooth[[i]]$fx)) {
        all(GamObj$smooth[[i]]$fx)
    } else {
        GamObj$smooth[[i]]$fixed
    }
    if (!fx && GamObj$smooth[[i]]$null.space.dim == 0 && !is.null(GamObj$R)) {
        ## random effect or fully penalized term
        res <- if (RanEff) {
            mgcv:::reTest(GamObj, i)
        } else {
            NULL
        }
        ## Test the mth smooth for equality to zero (m is not the RE term) and
        ## accounting for all random effects in model Inverted Nychka interval
        ## statistics
    } else {
        if (Quasi) {
            rdf <- resi_df
        } else {
            rdf <- -1
        }
        res <- mgcv:::testStat(p, Xt, V, min(ncol(Xt), edf1i), type = 0,
            res.df = rdf)
    }
    return(list(res = res, edfi = edfi, edf1i = edf1i))
}

#' @title Some checks for binomRegMethModelChecks
#'
#' @description Check if inputs fit one anoter according to
#' there shapes
#' @param data a data frame with rows as individual CpGs
#' appeared in all the samples.
#' The first 4 columns should contain the information of
#' `Meth_Counts` (methylated counts),
#' `Total_Counts` (read depths), `Position` (Genomic
#' position for the CpG site) and `ID`
#' (sample ID). The covariate information, such as
#' disease status or cell type composition
#' are listed in column 5 and onwards.
#' @param Z matrix for the covariate information
#' @param n.k a vector of basis dimensions for the
#' intercept and individual covariates.
#' \code{n.k} specifies an upper limit of the degrees
#' of each functional parameters.
#' @author Simon Laurin-Lemay, Kaiqiong Zhao
#' @noRd
binomRegMethModelChecks <- function(data, Z, n.k) {
    if (length(n.k) != (ncol(Z) + 1)) {
        stop("The length of n.k should equal to the number of covariates plus 1
             (for the intercept)")
    }
    if (any(data$X == 0)) {
        stop("The rows with Total_Counts equal to 0 should be deleted beforehand")
    }
    if (any(is.na(Z))) {
        stop("The covariate information should not have missing values")
    }
    if (any(!is.numeric(Z))) {
        stop("Please transform the covariate information into numeric values, eg.
             use dummy variables for the categorical covariates")
    }
}

#' @title Init for binomRegMethModel
#'
#' @description Initialize data
#' @param data a data frame with rows as individual CpGs
#' appeared in all the samples.
#' The first 4 columns should contain the information of
#' `Meth_Counts` (methylated counts),
#' `Total_Counts` (read depths), `Position` (Genomic position
#' for the CpG site) and `ID`
#' (sample ID). The covariate information, such as disease
#' status or cell type composition
#' are listed in column 5 and onwards.
#' @param covs  vector of covariate names. The covariates
#' with names in \code{covs} will
#' be included in the model and their covariate
#' effects will be estimated. The default is to fit all
#' covariates in \code{data}
#' @param n.k a vector of basis dimensions for the intercept
#' and individual covariates.
#' \code{n.k} specifies an upper limit of the degrees of each
#' functional parameters.
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{data}: the first four columns are 'Y', 'X',
#' 'Posit' and 'ID'; and the rest
#' are covariate values
#' \item \code{Z}: the matrix obtained by removing the first
#' 4 columns of data
#' }
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
binomRegMethModelInit <- function(data, covs, n.k) {
    data <- data.frame(data)
    if (is.factor(data$Position)) {
        ## message('The Position in the data set should be numeric other than a
        ## factor')
        data$Position <- as.numeric(as.character(data$Position))
    }
    if (any(!c("Meth_Counts", "Total_Counts", "Position") %in% colnames(data))) {
        stop("Please make sure object \"data\" have columns named as \"Meth_Counts\",
             \"Total_Counts\" and \"Position\" ")
    }
    colnames(data)[match(c("Meth_Counts", "Total_Counts", "Position"),
        colnames(data))] <- c("Y", "X", "Posit")
    if (is.null(covs)) {
        Z <- as.matrix(data[, -seq_len(4)], ncol = ncol(data) - 4)
        colnames(Z) <- colnames(data)[-seq_len(4)]
        binomRegMethModelChecks(data = data, Z = Z, n.k = n.k)
        return(out <- list(data = data, Z = Z))
    } else {
        id <- match(covs, colnames(data))
        if (any(is.na(id))) {
            stop(paste(covs[is.na(id)], " is not found in the input data frame"))
        } else {
            Z <- as.matrix(data[, id], ncol = length(id))
            colnames(Z) <- covs
            binomRegMethModelChecks(data = data, Z = Z, n.k = n.k)
            return(out <- list(data = data, Z = Z))
        }
    }
}

#' @title Fit gam for the initial value
#'
#' @description Fit gam for the initial value
#' @param data a data frame with rows as individual CpGs
#' appeared in all the samples.
#' The first 4 columns should contain the information of
#' `Meth_Counts` (methylated counts),
#' `Total_Counts` (read depths), `Position` (Genomic
#' position for the CpG site) and `ID`
#' (sample ID). The covariate information, such as
#' disease status or cell type composition
#' are listed in column 5 and onwards.
#' @param Quasi whether a Quasi-likelihood estimation
#' approach will be used
#' @param binom.link the link function used in the
#' binomial regression model; the
#' default is the logit link
#' @param method the method used to estimate the
#' smoothing parameters. The default
#' is the 'REML' method which is generally better than
#' prediction based criterion
#' \code{GCV.cp}
#' @param RanEff whether sample-level random effects
#' are added or not
#' @param scale nagative values mean scale paramter
#' should be estimated; if a positive
#' value is provided, a fixed scale will be used.
#' @param Z covariate matrix
#' @param n.k  a vector of basis dimensions for the
#' intercept and individual covariates.
#' \code{n.k} specifies an upper limit of the
#' degrees of each functional parameters
#' @return This function return a \code{list}
#' including objects:
#' \itemize{
#' \item \code{data}: a data frame with rows as individual
#' CpGs appeared
#' in all the
#' samples. The first 4 columns should contain the information
#' of `Meth_Counts`
#' (methylated counts), `Total_Counts` (read depths), `Position
#' (Genomic position
#' for the CpG site) and `ID` (sample ID). The covariate information,
#' such as disease
#' status or cell type composition are listed in column 5 and onwards.
#' \item \code{gam.int}: a gam object from mgcv::gam
#' \item \code{my.covar.fm}: the formula used to fit the gam
#' }
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
fitGam <- function(data, Quasi, binom.link, method, RanEff, scale, Z, n.k) {
    ## The smoothing formula corresponding to the Z
    formula.z.part <- vapply(seq_len(ncol(Z)), function(i) {
        paste0("s(Posit, k=n.k[", i + 1, "], fx=FALSE, bs=\"cr\", by=Z[,",
            i, "])")
    }, "")
    my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=FALSE, bs=\"cr\")", formula.z.part),
        collapse = "+")

    if (RanEff) {
        my.covar.fm <- paste0(my.covar.fm, "+ s(ID, bs=\"re\")")
        data$ID <- as.factor(data$ID)
    }
    ## Fit gam for the initial value
    if (Quasi) {
        gam.int <- mgcv::gam(as.formula(paste0("Y/X ~", my.covar.fm)),
            family = quasibinomial(link = binom.link), weights = data$X,
            data = data, method = method, scale = scale)
        return(out <- list(data = data, gam.int = gam.int, my.covar.fm = my.covar.fm))
    } else {
        gam.int <- mgcv::gam(as.formula(paste0("Y/X ~", my.covar.fm)),
            family = binomial(link = binom.link), weights = data$X, data = data,
            method = method, scale = scale)
        return(out <- list(data = data, gam.int = gam.int, my.covar.fm = my.covar.fm))
    }
}

#' @title phiFletcher calculate the Fletcher-based dispersion
#' estimate from the final gam output
#'
#' @description phiFletcher calculate the Fletcher-based dispersion
#' estimate from the
#' final gam output
#' @param data a data frame with rows as individual CpGs appeared
#' in all the samples.
#' The first 4 columns should contain the information of `Meth_Counts`
#' (methylated counts),
#' `Total_Counts` (read depths), `Position` (Genomic position for the
#' CpG site) and
#' `ID` (sample ID). The covariate information, such as disease
#' status or cell type
#' composition are listed in column 5 and onwards.
#' @param Quasi whether a Quasi-likelihood estimation approach
#' will be used
#' @param reml.scale whether a REML-based scale (dispersion)
#' estimator is used. The default
#' is Fletcher-based estimator
#' @param scale nagative values mean scale paramter should be estimated;
#' if a positive value
#' is provided, a fixed scale will be used.
#' @param gam.int a gam object from mgcv::gam
#' @return This function return a \code{phi_fletcher} object:
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
phiFletcher <- function(data, Quasi, reml.scale, scale, gam.int) {
    old.pi.ij <- gam.int$fitted.values
    edf.out <- gam.int$edf
    pearsonResiduals <- residuals(gam.int, type = "pearson")
    if (Quasi & scale <= 0) {
        if (reml.scale) {
            return(phi_fletcher <- gam.int$reml.scale)
        }
        ## * sqrt(data$X)
        my_s <- (1 - 2 * old.pi.ij)/(data$X * old.pi.ij * (1 - old.pi.ij)) *
            (data$Y - data$X * old.pi.ij)
        ## Note: the estimator implemented in the mgcv calculated my_s with an
        ## additional multiplier sqrt(data$X) But from the paper there shouldn't
        ## be this one phi_p=sum( pearsonResiduals^2)/(length(data$Y) -
        ## sum(edf.out))

        ## Feb 21, 2020 use edf1 to calculate the pearson's dispersion estimate
        ## not the edf; because I will use the asympototic chi-square dist. of
        ## phi.est March 3, 2020, the dispersion parameter estimated in gam use
        ## the edf.out instead of edf1.out
        phi_p <- sum(pearsonResiduals^2)/(length(data$Y) - sum(edf.out))
        return(phi_fletcher <- phi_p/(1 + mean(my_s)))
    }
    if (!Quasi) {
        return(phi_fletcher <- 1)
    }
    if (scale > 0) {
        return(phi_fletcher <- scale)
    }
}


#' @title extract design matrix
#'
#' @description extract design matrix
#' @param GamObj a mgcv object
#' @return This function return a design matrix \code{X_d}: R matrix
#' from the QR decomposition
#' @importFrom mgcv predict.gam
#' @importFrom mgcv model.matrix.gam
#' @author Kaiqiong Zhao Simon Laurin-Lemay
#' @noRd
extractDesignMatrix <- function(GamObj) {
    ## A more efficient way to extract design matrix. use a random sample of
    ## rows of the data to reduce the computational cost
    if (!is.null(GamObj$R)) {
        return(X_d <- GamObj$R)
    } else {
        sub.samp <- max(1000, 2 * length(GamObj$coefficients))
        if (nrow(GamObj$model) > sub.samp) {
            ## subsample to get X for p-values calc.  sample these rows from X
            ind <- sample(seq_len(nrow(GamObj$model)), sub.samp, replace = FALSE)
            X_d <- mgcv::predict.gam(GamObj, GamObj$model[ind, ], type = "lpmatrix")
        } else {
            ## don't need to subsample
            X_d <- mgcv::model.matrix.gam(GamObj)
        }
        ## exclude NA's (possible under na.exclude)
        return(X_d <- X_d[!is.na(rowSums(X_d)), ])
    }
}
