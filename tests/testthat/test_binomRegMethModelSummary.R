context("testing binomRegMethModelSummary")

### creation of the BEM.obj and BEM.obj.RA 
# data(RAdat)
# RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# BEM.obj <- binomRegMethModel(
#   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#   Quasi = FALSE, RanEff = FALSE, verbose = FALSE)
# BEM.obj.RA <- binomRegMethModel(
#   data=RAdat.f, n.k=rep(5, 2), p0=0.003307034, p1=0.9,
#   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#   Quasi = FALSE, RanEff = FALSE, covs = c("RA"), verbose = FALSE)

### creation of the inputs for binomRegMethModelSummary
### this code is coming from the binomRegMethModel function 
### just before calling the binomRegMethModelSummary function
# data=RAdat.f
# n.k=rep(5, 3)
# p0=0.003307034
# p1=0.9
# epsilon=10^(-6)
# epsilon.lambda=10^(-3)
# maxStep=200
# Quasi = FALSE
# RanEff = FALSE
# verbose = FALSE
# binom.link = "logit"
# method = "REML"
# covs = NULL
# reml.scale = FALSE
# scale = -2
# ### FOR BEM.obj
# initOut <- binomRegMethModelInit(data=data, covs=covs, n.k=n.k, 
#                                  verbose = verbose)
# Z <- initOut$Z
# fitGamOut <- fitGam(data=initOut$data, Quasi=Quasi, binom.link=binom.link,
#                     method=method, RanEff=RanEff, scale=scale, Z=Z, n.k=n.k,
#                     verbose = verbose)
# phi_fletcher <- phiFletcher(data=initOut$data, Quasi=Quasi, 
#                             reml.scale=reml.scale, scale=scale, 
#                             gam.int=fitGamOut$gam.int, verbose=verbose)
# fitEMOut <- fitEM(p0=p0, p1=p1, fitGamOut=fitGamOut, n.k=n.k, 
#                   binom.link=binom.link, method=method, Z=Z, Quasi=Quasi, 
#                   scale=scale, reml.scale=reml.scale, 
#                   phi_fletcher=phi_fletcher, epsilon=epsilon, 
#                   maxStep=maxStep, data=initOut$data,
#                   verbose = verbose)
# out <- fitEMOut$out
# Est.points <- fitEMOut$Est.points
# phi_fletcher <- out$phi_fletcher
# my.design.matrix <- mgcv::model.matrix.gam(out$GamObj)
# lengthUniqueDataID <- length(unique(initOut$data$ID))
# phi_reml <- out$GamObj$reml.scale
# estimateBZOut <- estimateBZ(Posit=initOut$data$Posit, 
#                             my.design.matrix=my.design.matrix, 
#                             ncolsZ=ncol(Z), n.k=n.k, verbose = verbose)
# Beta.out <- estimateBeta(BZ=estimateBZOut$BZ, BZ.beta=estimateBZOut$BZ.beta,
#                          n.k=n.k, Z=Z, out=out, verbose = verbose)
# estimateVarOut <- estimateVar(out=out, phi_fletcher=phi_fletcher, 
#                               data=initOut$data, 
#                               my.design.matrix=my.design.matrix, Z=Z, p0=p0,
#                               p1=p1, RanEff=RanEff, 
#                               lengthUniqueDataID=lengthUniqueDataID,n.k=n.k,
#                               verbose = verbose)
# estimateSEOut <- estimateSE(estimateBZOut=estimateBZOut, Z=Z, 
#                             estimateVarOut=estimateVarOut, 
#                             phi_fletcher=phi_fletcher, phi_reml=phi_reml,
#                             verbose = verbose)
# X_d <- extractDesignMatrix(GamObj=out$GamObj, verbose = verbose)
# ## and p value calculation tr(2A - A^2)
# edf1.out <- out$edf1
# ## Effective degree of freedom: edf --trace of the hat matrix
# edf.out <- out$edf
# resi_df <- nrow(initOut$data) - sum(edf.out)
# 
# inputs <- list(GamObj=out$GamObj, 
#                var.cov.alpha=estimateVarOut$var.cov.alpha,
#                new.par=out$par, edf.out=edf.out, 
#                edf1.out=edf1.out, X_d=X_d, 
#                resi_df=resi_df, Quasi=Quasi, 
#                scale=scale, RanEff=RanEff, Z=Z,
#                verbose = verbose)
# ### FOR BEM.obj.RA
# covs <- c("RA")
# n.k <- rep(5, 2)
# initOut <- binomRegMethModelInit(data=data, covs=covs, n.k=n.k,
#                                  verbose = verbose)
# Z <- initOut$Z
# fitGamOut <- fitGam(data=initOut$data, Quasi=Quasi, binom.link=binom.link,
#                     method=method, RanEff=RanEff, scale=scale, Z=Z, n.k=n.k,
#                     verbose = verbose)
# phi_fletcher <- phiFletcher(data=initOut$data, Quasi=Quasi,
#                             reml.scale=reml.scale, scale=scale,
#                             gam.int=fitGamOut$gam.int, verbose=verbose)
# fitEMOut <- fitEM(p0=p0, p1=p1, fitGamOut=fitGamOut, n.k=n.k,
#                   binom.link=binom.link, method=method, Z=Z, Quasi=Quasi,
#                   scale=scale, reml.scale=reml.scale,
#                   phi_fletcher=phi_fletcher, epsilon=epsilon,
#                   maxStep=maxStep, data=initOut$data,
#                   verbose = verbose)
# out <- fitEMOut$out
# Est.points <- fitEMOut$Est.points
# phi_fletcher <- out$phi_fletcher
# my.design.matrix <- mgcv::model.matrix.gam(out$GamObj)
# lengthUniqueDataID <- length(unique(initOut$data$ID))
# phi_reml <- out$GamObj$reml.scale
# estimateBZOut <- estimateBZ(Posit=initOut$data$Posit,
#                             my.design.matrix=my.design.matrix,
#                             ncolsZ=ncol(Z), n.k=n.k, verbose = verbose)
# Beta.out <- estimateBeta(BZ=estimateBZOut$BZ, BZ.beta=estimateBZOut$BZ.beta,
#                          n.k=n.k, Z=Z, out=out, verbose = verbose)
# estimateVarOut <- estimateVar(out=out, phi_fletcher=phi_fletcher,
#                               data=initOut$data,
#                               my.design.matrix=my.design.matrix, Z=Z, p0=p0,
#                               p1=p1, RanEff=RanEff,
#                               lengthUniqueDataID=lengthUniqueDataID,n.k=n.k,
#                               verbose = verbose)
# estimateSEOut <- estimateSE(estimateBZOut=estimateBZOut, Z=Z,
#                             estimateVarOut=estimateVarOut,
#                             phi_fletcher=phi_fletcher, phi_reml=phi_reml,
#                             verbose = verbose)
# X_d <- extractDesignMatrix(GamObj=out$GamObj, verbose = verbose)
# ## and p value calculation tr(2A - A^2)
# edf1.out <- out$edf1
# ## Effective degree of freedom: edf --trace of the hat matrix
# edf.out <- out$edf
# resi_df <- nrow(initOut$data) - sum(edf.out)
# 
# inputs.RA <- list(GamObj=out$GamObj,
#                var.cov.alpha=estimateVarOut$var.cov.alpha,
#                new.par=out$par, edf.out=edf.out,
#                edf1.out=edf1.out, X_d=X_d,
#                resi_df=resi_df, Quasi=Quasi,
#                scale=scale, RanEff=RanEff, Z=Z,
#                verbose = verbose)

load(paste0(getwd(), "/data/BEM.obj.rda"))
load(paste0(getwd(), "/data/BEM.obj.RA.rda"))
load(paste0(getwd(), "/data/ref_input_binomRegMethModelSummary.rda"))
load(paste0(getwd(), "/data/ref_input_binomRegMethModelSummary.RA.rda"))

s.table <- binomRegMethModelSummary(inputs$GamObj, inputs$var.cov.alpha,
                                    inputs$new.par ,inputs$edf.out,
                                    inputs$edf1.out, inputs$X_d,
                                    inputs$resi_df, inputs$Quasi,
                                    inputs$scale, inputs$RanEff,
                                    inputs$Z, verbose = FALSE)

s.table.RA <- binomRegMethModelSummary(inputs.RA$GamObj, inputs.RA$var.cov.alpha,
                                    inputs.RA$new.par ,inputs.RA$edf.out,
                                    inputs.RA$edf1.out, inputs.RA$X_d,
                                    inputs.RA$resi_df, inputs.RA$Quasi,
                                    inputs.RA$scale, inputs.RA$RanEff,
                                    inputs.RA$Z, verbose = FALSE)

test_that("binomRegMethModelSummary provides the correct regional p-value summary", {
  expect_true(all(c("Intercept","T_cell", "RA") %in% rownames(s.table)))
  expect_true(all(c("Intercept", "RA") %in% rownames(s.table.RA)))
  expect_false(all(c("Intercept","T_cell", "RA") %in% rownames(s.table.RA)))
  expect_true(all.equal(s.table, BEM.obj$reg.out, tolerance = 1e-3))
  expect_true(all.equal(s.table.RA, BEM.obj.RA$reg.out, tolerance = 1e-3))
})




