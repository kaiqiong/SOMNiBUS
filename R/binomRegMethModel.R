#' @title A smoothed-EM algorithm to estimate covariate effects and test regional association in Bisulfite Sequencing-derived methylation data
#'
#' @description This function fits a (dispersion-adjusted) binomial regression model to regional methylation data,
#' and reports the estimated smooth covariate effects and regional p-values for the test of DMRs
#' (differentially methylation regions). Over or under dispersion across loci is accounted for in the model
#' by the combination of a multiplicative dispersion parameter (or scale parameter) and a sample-specific random effect.
#' @description This method can deal with outcomes, i.e. the number of methylated reads in a region, that are contaminated
#'  by known false methylation calling rate (\code{p0}) and false non-methylation calling rate (\code{1-p1}).
#' @description The covariate effects are assumed to smoothly vary across genomic regions. In order to estimate them,
#' the algorithm first represents the functional paramters by a linear combination of a set of restricted cubic splines
#' (with dimention \code{n.k}), and a smoothness penalization term which depends on the smoothing parameters \code{lambdas}
#' is also added to control smoothness. The estimation is performed by an iterated EM algorithm. Each M step constitutes
#' an outer Newton's iteration to estimate smoothing parameters \code{lambdas} and an inner P-IRLS iteration to estimate spline
#' coefficients \code{alpha} for the covariate effects. Currently, the computation in the M step depends the implementation of
#' \code{gam()} in package \code{mgcv}.
#' @param data a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the
#' information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site)
#' and `ID` (sample ID). The covariate information, such as disease status or cell type composition, are listed in column 5 and onwards.
#' @param n.k a vector of basis dimensions for the intercept and individual covariates. \code{n.k} specifies an upper limit of the
#' degrees of each functional parameters.
#' @param p0 the probability of observing a methylated read when the underlying true status is unmethylated. \code{p0} is the rate
#' of false methylation calls, i.e. false positive rate.
#' @param p1 the probability of observing a methylated read when the underlying true status is methylated. \code{1-p1} is the rate
#' of false non-methylation calls, i.e. false negative rate.
#' @param covs a vector of covariate names. The covariates with names in \code{covs} will be included in the model and their covariate
#' effects will be estimated. The default is to fit all covariates in \code{data}
#' @param Quasi whether a Quasi-likelihood estimation approach will be used; in other words, whether a multiplicative dispersion is added
#' in the model or not.
#' @param epsilon numeric; stopping criterion for the closeness of estimates of spline coefficients from two consecutive iterations.
#' @param epsilon.lambda numeric; stopping criterion for the closeness of estimates of smoothing parameter \code{lambda} from two consecutive
#' iterations.
#' @param maxStep the algorithm will step if the iteration steps exceed \code{maxStep}
#' @param detail indicate whether print the number of iterations
#' @param binom.link the link function used in the binomial regression model; the default is the logit link
#' @param method the method used to estimate the smoothing parameters. The default is the 'REML' method which is generally better than
#' prediction based criterion \code{GCV.cp}
#' @param RanEff whether sample-level random effects are added or not
#' @param reml.scale whether a REML-based scale (dispersion) estimator is used. The default is Fletcher-based estimator
#' @param scale nagative values mean scale paramter should be estimated; if a positive value is provided, a fixed scale will be used.
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{est}: estimates of the spline basis coefficients \code{alpha}
#' \item \code{lambda}: estimates of the smoothing parameters for each functional paramters
#' \item \code{est.pi}: predicted methylation lelves for each row in the input \code{data}
#' \item \code{ite.points}: estimates of \code{est}, \code{lambda} at each EM iteration
#' \item \code{cov1}: estimated variance-covariance matrix of the basis coefficients alphs
#' \item \code{reg.out}: regional testing output obtained using Fletcher-based dispersion estimate;
#' an additional "ID" row would appear if RanEff is TRUE
#' \item \code{reg.out.reml.scale}:regional testing output obtained using REML-based dispersion estimate;
#' \item \code{reg.out.gam}:regional testing output obtained using (Fletcher-based) dispersion estimate from mgcv package;
#' \item \code{phi_fletcher}: Fletcher-based estimate of the (multiplicative) dispersion parameter
#' \item \code{phi_reml}: REML-based estimate of the (multiplicative) dispersion parameter
#' \item \code{SE.out}: a matrix of the estimated pointwise Standard Errors (SE); number of rows are the number of unique CpG
#' sites in the input data and the number of columns equal to the total number of covariates fitted in the model
#' (the first one is the intercept)
#' \item \code{SE.out.REML.scale}:
#' \item \code{uni.pos}: the genomic postions for each row of CpG sites in the matrix \code{SE.out}
#' \item \code{Beta.out}: a matrix of the estimated covariate effects beta(t), here t denots the genomic positions.
#' \item \code{ncovs}: number of functional paramters in the model (including the intercept)
#' \item \code{sigma00}: estimated variance for the random effect if RanEff is TRUE; NA if RanEff is FALSE
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
#'   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200, detail=FALSE
#' )
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom mgcv model.matrix.gam
#' @importFrom mgcv s
#' @importFrom Matrix bdiag
#' @importFrom stats as.formula binomial pchisq rbinom rnorm quasibinomial residuals predict model.matrix runif
#' @export
binomRegMethModel <- function(data, n.k, p0=0.003, p1=0.9, Quasi=TRUE, epsilon=10^(-6),
    epsilon.lambda=10^(-3), maxStep=200, detail=FALSE, binom.link="logit",
    method="REML", covs=NULL, RanEff=TRUE, reml.scale=FALSE, scale=-2) {
    ## an error of 'n.k is not found' would appear if without this global
    ## environment assignment; so I save n.k in a parent scope
    n.k <<- n.k

    initOut<-binomRegMethModelInit(data=data, covs=covs, n.k=n.k)
    Z<-initOut$Z
    fitGamOut<-fitGam(data=initOut$data, Quasi=Quasi, binom.link=binom.link, method=method, RanEff=RanEff, scale=scale, Z=Z)


    ## Estimates
    phi_fletcher <- phiFletcher(data=fitGamOut$data, Quasi=Quasi, reml.scale=reml.scale, scale=scale, gam.int=fitGamOut$gam.int)
    if (p0 > 0 | p1 < 1) {
        ## code used to generate
        ## /tests/testthat/data/ref_input_binomRegMethModelUpdate.RDS input =
        ## list(data=data, old.pi.ij=old.pi.ij, p0=p0, p1=p1, n.k=n.k,
        ## binom.link=binom.link, method=method, Z=Z, my.covar.fm=my.covar.fm,
        ## Quasi=Quasi, scale=phi_fletcher) path_ref_input_binomRegMethModelUpdate <-
        ## paste(paste(getwd(), '/tests/testthat/data/', sep=''),
        ## 'ref_input_binomRegMethModelUpdate.RDS', sep='') saveRDS(input,
        ## path_ref_input_binomRegMethModelUpdate)

        out <- binomRegMethModelUpdate(data=fitGamOut$data, pi.ij=fitGamOut$gam.int$fitted.values, p0=p0, p1=p1, n.k=n.k,
            binom.link=binom.link, method=method, Z=Z, my.covar.fm=fitGamOut$my.covar.fm,
            Quasi=Quasi, scale=phi_fletcher)
        Est.points <- rbind(c(fitGamOut$gam.int$coefficients, fitGamOut$gam.int$sp, phi_fletcher), c(out$par, out$lambda, out$phi_fletcher))
        ## Do the iteration The stopping criterion is that estimator a and b are
        ## close enough I exclude the criterio that lambda-a and lambda-b are
        ## close
        old.pi.ij <- fitGamOut$gam.int$fitted.values
        iter <- 1
        while (sqrt(sum((out$pi.ij - old.pi.ij)^2)) > epsilon & iter < maxStep) {
            if (detail) {
                message(paste0("iteration", iter))
            }
            iter <- iter + 1
            old.pi.ij <- out$pi.ij
            out <- binomRegMethModelUpdate(data=fitGamOut$data, pi.ij=old.pi.ij, p0=p0, p1=p1, n.k=n.k, binom.link=binom.link,
                method=method, Z=Z, my.covar.fm=fitGamOut$my.covar.fm, Quasi=Quasi,
                scale=phi_fletcher)
            Est.points <- rbind(Est.points, c(out$par,out$lambda, out$phi_fletcher))
        }
    } else {
        out <- list(pi.ij=fitGamOut$gam.int$fitted.values, par=fitGamOut$gam.int$coefficients,
                    lambda=fitGamOut$gam.int$sp, edf1=fitGamOut$gam.int$edf1,
                    pearson_res=residuals(fitGamOut$gam.int, type="pearson"),
                    deviance_res=residuals(fitGamOut$gam.int, type ="deviance"),
                    edf=fitGamOut$gam.int$edf, phi_fletcher=phi_fletcher, GamObj=fitGamOut$gam.int)
        Est.points <- c(fitGamOut$gam.int$coefficients, fitGamOut$gam.int$sp, phi_fletcher)
    }

    ## Estimated dispersion paramters (Fletcher adjusted)
    phi_fletcher <- out$phi_fletcher

    ##-------------------------------------------
    ## Calculate SE of alpha
    ##------------------------------------------
    ## the part to calculate the SE of alpha the model matrix for the GamObj
    ## and the FinalGamObj is the same. the difference is only on the
    ## outcomes
    my.design.matrix <- mgcv::model.matrix.gam(out$GamObj)
    lengthUniqueDataID <- length(unique(fitGamOut$data$ID))
    phi_reml <- out$GamObj$reml.scale
    ##------------------------------------------------------------------------
    ## from the variance-covariance of alpha, to report 1. var(beta_p(t)) 2.
    ## Report the chi-square statistic for each covariate 3. Report the
    ## p-value for each covariate
    ##------------------------------------------------------------------------

    estimateBZOut <- estimateBZ(data=fitGamOut$data, my.design.matrix=my.design.matrix, ncolsZ = ncol(Z), n.k=n.k)
    Beta.out <- estimateBeta(BZ=estimateBZOut$BZ, BZ.beta=estimateBZOut$BZ.beta, n.k=n.k, Z=Z, out=out)
    cum_s <- cumsum(n.k)

    ##---------------------------------------------------------------
    ## calculate var_cov (for alpha & beta) from the hessianComp matrix
    ##---------------------------------------------------------------
    w_ij <- out$pi.ij * (1 - out$pi.ij)/phi_fletcher
    H <- hessianComp(w_ij=w_ij, new.par=out$par,
        new.lambda=out$lambda,  X=fitGamOut$data$X, Y=fitGamOut$data$Y, my.design.matrix=my.design.matrix, gam.int=out$GamObj, Z=Z, pred.pi=out$pi.ij, p0=p0, p1=p1,
        disp_est=phi_fletcher, RanEff=RanEff, N=lengthUniqueDataID)

    var.cov.alpha <- solve(-H)
    var.alpha.0 <- var.cov.alpha[seq_len(n.k[1]), seq_len(n.k[1])]
    var.alpha.sep <- lapply(seq_len(ncol(Z)), function(i) {
        var.cov.alpha[(cum_s[i] + 1):cum_s[i + 1], (cum_s[i] + 1):cum_s[i +
            1]]
    })

    SE.out <- cbind(sqrt(pmax(0, rowSums((estimateBZOut$BZ %*% var.alpha.0) * estimateBZOut$BZ))),
        vapply(seq_len(ncol(Z)), function(i) {
            sqrt(pmax(0, rowSums((estimateBZOut$BZ.beta[[i]] %*% var.alpha.sep[[i]]) *
                estimateBZOut$BZ.beta[[i]])))
        }, FUN.VALUE = rep(1.0, length(estimateBZOut$uni.pos))))

    rownames(SE.out) <- estimateBZOut$uni.pos
    colnames(SE.out) <- c("Intercept", colnames(Z))
    SE.pos <- estimateBZOut$uni.pos
    SE.out.REML.scael <- SE.out/sqrt(phi_fletcher) * sqrt(phi_reml)

    ##---------------------------------------------------------------
    ## calculate the region-based statistic from the testStats function in
    ## mgcv
    ##---------------------------------------------------------------

    ## p : estimated paramters --- alpha.0, alpha.seq Xt: the design matrix
    ## --- my.design.matrix V: estimated variance matrix on entry `rank'
    ## should be an edf estimate 0. Default using the fractionally truncated
    ## pinv.  1. Round down to k if k<= rank < k+0.05, otherwise up.  res.df
    ## is residual dof used to estimate scale. <=0 implies fixed scale.

    X_d<-extractDesignMatrix(GamObj=out$GamObj)

    ## Effective degrees of freedom: edf1 -- good for chisquare test and p
    ## value calculation tr(2A - A^2)
    edf1.out <- out$edf1
    ## Effective degree of freedom: edf --trace of the hat matrix
    edf.out <- out$edf
    ## the residuals degrees of freedom
    resi_df <- nrow(fitGamOut$data) - sum(edf.out)
    s.table <- binomRegMethModelSummary(out$GamObj, var.cov.alpha, out$par, edf.out,
        edf1.out, X_d, resi_df, Quasi, scale, RanEff, Z)
    s.table.REML.scale <- binomRegMethModelSummary(out$GamObj, var.cov.alpha/phi_fletcher *
        phi_reml, out$par, edf.out, edf1.out, X_d, resi_df, Quasi, scale,
        RanEff, Z)

    ## var_out=list(cov1=var.cov.alpha, reg.out=reg.out, SE.out=SE.out,
    ## uni.pos=SE.pos, pvalue=pvalue , ncovs=ncol(Z)+1) Est_out=list(est =
    ## new.par, lambda=new.lambda, est.pi=new.pi.ij, ite.points=Est.points,
    ## Beta.out=Beta.out, phi_fletcher=phi_fletcher)


    if (RanEff) {
        sigma00 <- out$GamObj$reml.scale/out$GamObj$sp["s(ID)"]
    } else {
        sigma00 <- NA
    }

    reg.out.gam <- summary(out$GamObj)$s.table

    return(out<-list(est=out$par, lambda=out$lambda, est.pi=out$pi.ij,
        Beta.out=Beta.out, phi_fletcher=phi_fletcher, phi_reml=phi_reml,
        reg.out=s.table, reg.out.reml.scale=s.table.REML.scale, cov1=var.cov.alpha,
        reg.out.gam=reg.out.gam, SE.out=SE.out, SE.out.REML.scale=SE.out.REML.scael,
        uni.pos=SE.pos, ncovs=ncol(Z) + 1, ite.points=Est.points,
        sigma00=sigma00))
}
#' @title Get the basis matrix for beta0(t) and beta(t)
#'
#' @description This function aims to extract the basis matrix BZ, and BZ.beta from the expansion beta0(t) = BZ * alpha0
#' and beta(t) = BZ.beta * alpha.
#' @param data a data frame with one column named "Posit"
#' @param my.design.matrix the design matrix for a \code{fitGamOut} with data as input and Z as covariates.
#' @param ncolsZ number of columns of covariate matrix Z

#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{uni.pos} the unique genomic positions present in the data
#' \item \code{BZ} basis matrix for the intercept beta0(t); a matrix with \code{length(uni.pos)} rows and \code{length(alpha0)} columns
#' \item \code{BZ.beta} a list of length \code{ncolsZ}; each corresponds to the basis matrix for beta_p(t), p = 1, 2, .. \code{ncolsZ}
#' }
#' @author  Kaiqiong Zhao, Simon Laurin-Lemay
#' @import mgcv
#' @noRd
estimateBZ <- function(data, my.design.matrix, ncolsZ, n.k){
    uni.pos <- unique(data$Posit)
    uni.id <- match(uni.pos, data$Posit)
    BZ <- my.design.matrix[uni.id, seq_len(n.k[1])]
    BZ.beta <- lapply(seq_len(ncolsZ), function(i) {
        mgcv::smooth.construct(mgcv::s(Posit, k=n.k[i + 1], fx=FALSE,
            bs="cr"), data=data[uni.id, ], knots=NULL)$X
    })
    return(out<-list(uni.pos=uni.pos, BZ=BZ, BZ.beta=BZ.beta))
}

#' @title estimate of beta(t)
#'
#' @description Given a final GAM output, extract the estimates of beta(t) = BZ * alpha, where t is the unique genomic positions
#' present in the input data; only extract the fixed effect estimates
#' @param BZ basis matrix for the intercept beta0(t); a matrix with \code{length(uni.pos)} rows and \code{length(alpha0)} columns
#' @param BZ.beta a list of length \code{ncolsZ}; each corresponds to the basis matrix for beta_p(t)
#' @param n.k a vector of basis dimensions for the intercept and individual covariates. \code{n.k} specifies an upper limit of the degrees of each functional parameters.
#' @param Z a covariate matrix
#' @param out an object from \code{binomRegMethModelUpdate}
#' @return \code{Beta.out}: a matrix of the estimated covariate effects beta(t), here t denots the unique genomic positions for intercept
#' and all covariates. Beta.out = cbind(beta.0(t), beta.1(t), ... beta.P(t)); \code{length(uni.pos)} rows and \code{ncol(Z)+1} columns
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
estimateBeta <- function(BZ, BZ.beta, n.k, Z, out){
    ## ------ 3: estimate of beta(t) --- #
    cum_s <- cumsum(n.k)
    alpha.sep <- lapply(seq_len(ncol(Z)), function(i) {
        out$par[(cum_s[i] + 1):cum_s[i + 1]]
    })
    alpha.0 <- out$par[seq_len(n.k[1])]

    Beta.out <- cbind(BZ %*% alpha.0, vapply(seq_len(ncol(Z)), function(i) {
        BZ.beta[[i]] %*% alpha.sep[[i]]
    },FUN.VALUE = rep(1.0, length(estimateBZOut$uni.pos))))

    colnames(Beta.out) <- c("Intercept", colnames(Z))
    return (Beta.out)
}

#' @title binomRegMethModelSummary Extract the regional testing results from a GamObj
#'
#' @description Extract the regional testing results from a GamObj
#' @param GamObj an fit from mgcv::gam
#' @param var.cov.alpha Estimated variance-covariance matrix for coefficients alpha
#' @param new.par Estimated coefficients alpha
#' @param edf.out Effective degrees of freedom for each alpha, calculated from the trace of the hat matrix
#' @param edf1.out Effective degrees of freedom 2, calculated from tr(2A - A^2); good for chisquare test
#' @param X_d R from the QR decomposition for the design matrix
#' @param resi_df Residual degrees of freedom; nrow(data) - sum(edf.out)
#' @param Quasi Whether a multiplicative dispersion parameter is added or not; quasibinomial or binomial
#' @param scale nagative values mean scale paramter should be estimated; if a positive value is provided, a fixed scale will be used.
#' @param RanEff Whether a subject random effect is added or not
#' @param Z covariate matrix; \code{Z = data[,-c(1:4)]}
#' @return s.table a table for each covariate
#' @author  Kaiqiong Zhao
#' @import mgcv
#' @noRd
binomRegMethModelSummary <- function(GamObj, var.cov.alpha, new.par, edf.out, edf1.out,
    X_d, resi_df, Quasi, scale, RanEff, Z) {
    ii <- 0
    m <- length(GamObj$smooth)
    df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
    for (i in seq_len(m)) {
        start <- GamObj$smooth[[i]]$first.para
        stop <- GamObj$smooth[[i]]$last.para

        V <- var.cov.alpha[start:stop, start:stop, drop=FALSE]  ## Bayesian

        p <- new.par[start:stop]  ## params for smooth
        edfi <- sum(edf.out[start:stop])  ## edf for this smooth
        ## extract alternative edf estimate for this smooth, if possible...
        edf1i <- sum(edf1.out[start:stop])
        Xt <- X_d[, start:stop, drop=FALSE]
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
            res <- mgcv:::testStat(p, Xt, V, min(ncol(Xt), edf1i), type=0,
                res.df=rdf)
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

#' @title Some checks for binomRegMethModelChecks
#'
#' @description Check if inputs fit one anoter according to there shapes
#' @param data a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site) and `ID` (sample ID). The covariate information, such as disease status or cell type composition are listed in column 5 and onwards.
#' @param Z Lorem ipsum dolor sit amet
#' @param n.k a vector of basis dimensions for the intercept and individual covariates. \code{n.k} specifies an upper limit of the degrees of each functional parameters.
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
binomRegMethModelChecks <- function(data, Z, n.k){
    if (length(n.k) != (ncol(Z) + 1)) {
        stop("The length of n.k should equal to the number of covariates plus 1 (for the intercept)")
    }
    if (any(data$X == 0)) {
        stop("The rows with Total_Counts equal to 0 should be deleted beforehand")
    }
    if (any(is.na(Z))) {
        stop("The covariate information should not have missing values")
    }
    if (any(!is.numeric(Z))) {
        stop("Please transform the covariate information into numeric values, eg. use dummy variables for the categorical covariates")
    }
}

#' @title Init for binomRegMethModel
#'
#' @description Initialize data
#' @param data a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site) and `ID` (sample ID). The covariate information, such as disease status or cell type composition are listed in column 5 and onwards.
#' @param covs  vector of covariate names. The covariates with names in \code{covs} will be included in the model and their covariate
#' effects will be estimated. The default is to fit all covariates in \code{data}
#' @param n.k a vector of basis dimensions for the intercept and individual covariates. \code{n.k} specifies an upper limit of the degrees of each functional parameters.
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{data}: the first four columns are "Y", "X", "Posit" and "ID"; and the rest are covariate values
#' \item \code{Z}: the matrix obtained by removing the first 4 columns of data
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
        stop("Please make sure object \"data\" have columns named as \"Meth_Counts\", \"Total_Counts\" and \"Position\" ")
    }
    colnames(data)[match(c("Meth_Counts", "Total_Counts", "Position"),
        colnames(data))] <- c("Y", "X", "Posit")
    if (is.null(covs)) {
        Z <- as.matrix(data[, -seq_len(4)], ncol=ncol(data) - 4)
        colnames(Z) <- colnames(data)[-seq_len(4)]
        binomRegMethModelChecks(data=data, Z=Z, n.k=n.k)
        return(out<-list(data=data, Z=Z))
    } else {
        id <- match(covs, colnames(data))
        if (any(is.na(id))) {
            stop(paste(covs[is.na(id)], " is not found in the input data frame"))
        } else {
            Z <- as.matrix(data[, id], ncol=length(id))
            colnames(Z) <- covs
            binomRegMethModelChecks(data=data, Z=Z, n.k=n.k)
            return(out<-list(data=data, Z=Z))
        }
    }
}

#' @title Fit gam for the initial value
#'
#' @description Fit gam for the initial value
#' @param data a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site) and `ID` (sample ID). The covariate information, such as disease status or cell type composition are listed in column 5 and onwards.
#' @param Quasi whether a Quasi-likelihood estimation approach will be used
#' @param binom.link the link function used in the binomial regression model; the default is the logit link
#' @param method the method used to estimate the smoothing parameters. The default is the 'REML' method which is generally better than prediction based criterion \code{GCV.cp}
#' @param RanEff whether sample-level random effects are added or not
#' @param scale nagative values mean scale paramter should be estimated; if a positive value is provided, a fixed scale will be used.
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{data}: a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site) and `ID` (sample ID). The covariate information, such as disease status or cell type composition are listed in column 5 and onwards.
#' \item \code{gam.int}: a gam object from mgcv::gam
#' \item \code{my.covar.fm}: the formula used to fit the gam
#' }
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
fitGam<-function(data, Quasi, binom.link, method, RanEff, scale, Z){
    ## The smoothing formula corresponding to the Z
    formula.z.part <- vapply(seq_len(ncol(Z)), function(i) {
        paste0("s(Posit, k=n.k[", i + 1, "], fx=FALSE, bs=\"cr\", by=Z[,",
            i, "])")
    }, "")
    my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=FALSE, bs=\"cr\")", formula.z.part),
        collapse="+")

    if (RanEff) {
        my.covar.fm <- paste0(my.covar.fm, "+ s(ID, bs=\"re\")")
        data$ID <- as.factor(data$ID)
    }
    ## Fit gam for the initial value
    if (Quasi) {
        gam.int <- mgcv::gam(as.formula(paste0("Y/X ~", my.covar.fm)),
            family=quasibinomial(link=binom.link), weights=data$X, data=data,
            method=method, scale=scale)
        return(out<-list(data=data, gam.int=gam.int, my.covar.fm=my.covar.fm))
    } else {
        gam.int <- mgcv::gam(as.formula(paste0("Y/X ~", my.covar.fm)),
            family=binomial(link=binom.link), weights=data$X, data=data,
            method=method, scale=scale)
        return(out<-list(data=data, gam.int=gam.int, my.covar.fm=my.covar.fm))
    }
}

#' @title phiFletcher calculate the Fletcher-based dispersion estimate from the final gam output
#'
#' @description phiFletcher calculate the Fletcher-based dispersion estimate from the final gam output
#' @param data a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site) and `ID` (sample ID). The covariate information, such as disease status or cell type composition are listed in column 5 and onwards.
#' @param Quasi whether a Quasi-likelihood estimation approach will be used
#' @param reml.scale whether a REML-based scale (dispersion) estimator is used. The default is Fletcher-based estimator
#' @param scale nagative values mean scale paramter should be estimated; if a positive value is provided, a fixed scale will be used.
#' @param gam.int a gam object from mgcv::gam
#' @return This function return a \code{phi_fletcher} object:
#' @author Kaiqiong Zhao, Simon Laurin-Lemay
#' @noRd
phiFletcher<-function(data, Quasi, reml.scale, scale, gam.int){
    old.pi.ij<-gam.int$fitted.values
    old.par<-gam.int$coefficients
    lambda<-gam.int$sp
    edf.out<-gam.int$edf
    edf1.out<-gam.int$edf1
    pearsonResiduals<-residuals(gam.int, type="pearson")
    devianceResiduals<-residuals(gam.int, type="deviance")
    if (Quasi & scale <= 0) {
        if (reml.scale) {
            return(phi_fletcher <- gam.int$reml.scale)
        }
        ## * sqrt(data$X)
        my_s <- (1 - 2 * old.pi.ij)/(data$X * old.pi.ij * (1 - old.pi.ij)) *
            (data$Y - data$X * old.pi.ij)
        ## Note: the estimator implemented in the mgcv calculated my_s with an
        ## additional multiplier sqrt(data$X) But from the paper there shouldn't
        ## be this one phi_p=sum( pearsonResiduals^2)/(length(data$Y) - sum(edf.out))

        ## Feb 21, 2020 use edf1 to calculate the pearson's dispersion estimate
        ## not the edf; because I will use the asympototic chi-square dist. of
        ## phi.est March 3, 2020, the dispersion parameter estimated in gam use
        ## the edf.out instead of edf1.out
        phi_p <- sum(pearsonResiduals^2)/(length(data$Y) - sum(edf.out))
        return(phi_fletcher <- phi_p/(1 + mean(my_s)))
    }
    if (!Quasi) {
        return(phi_fletcher<-1)
    }
    if (scale > 0) {
        return(phi_fletcher<-scale)
    }
}


#' @title extract design matrix
#'
#' @description extract design matrix
#' @param GamObj a mgcv object
#' @return This function return a design matrix \code{X_d}: R matrix from the QR decomposition
#' @author Kaiqiong Zhao Simon Laurin-Lemay
#' @noRd
extractDesignMatrix<-function(GamObj){
    ## A more efficient way to extract design matrix. use a random sample of
    ## rows of the data to reduce the computational cost
        if (!is.null(GamObj$R)) {
            X_d <- GamObj$R
        } else {
            sub.samp <- max(1000, 2 * length(GamObj$coefficients))
            if (nrow(GamObj$model) > sub.samp) {
                ## subsample to get X for p-values calc.  sample these rows from X
                ind <- sample(seq_len(nrow(GamObj$model)), sub.samp, replace=FALSE)
                X_d <- predict(GamObj, GamObj$model[ind, ], type="lpmatrix")
            } else {
                ## don't need to subsample
                X_d <- model.matrix(GamObj)
            }
            ## exclude NA's (possible under na.exclude)
            return(X_d <- X_d[!is.na(rowSums(X_d)), ])
        }  ## end if (m>0)
    }
