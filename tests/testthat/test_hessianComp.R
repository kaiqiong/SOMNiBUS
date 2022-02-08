context("testing hessianComp")

### creation of the BEM.obj
data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# BEM.obj <- binomRegMethModel(
#   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#   Quasi = FALSE, RanEff = FALSE, verbose = FALSE)
load(paste0(getwd(), "/data/BEM.obj.rda"))

### creation of the inputs for hessianComp
### this code is coming from the binomRegMethModel and estimateVar functions 
### just before calling the estimateVar function
data=RAdat.f
n.k=rep(5, 3)
p0=0.003307034
p1=0.9
epsilon=10^(-6)
epsilon.lambda=10^(-3)
maxStep=200
Quasi = FALSE
RanEff = FALSE
verbose = FALSE
binom.link = "logit"
method = "REML"
covs = NULL
reml.scale = FALSE
scale = -2
initOut <- binomRegMethModelInit(data=data, covs=covs, n.k=n.k,
                                 verbose = verbose)
Z <- initOut$Z
fitGamOut <- fitGam(data=initOut$data, Quasi=Quasi, binom.link=binom.link,
                    method=method, RanEff=RanEff, scale=scale, Z=Z, n.k=n.k,
                    verbose = verbose)
phi_fletcher <- phiFletcher(data=initOut$data, Quasi=Quasi,
                            reml.scale=reml.scale, scale=scale,
                            gam.int=fitGamOut$gam.int, verbose=verbose)
fitEMOut <- fitEM(p0=p0, p1=p1, fitGamOut=fitGamOut, n.k=n.k,
                  binom.link=binom.link, method=method, Z=Z, Quasi=Quasi,
                  scale=scale, reml.scale=reml.scale,
                  phi_fletcher=phi_fletcher, epsilon=epsilon,
                  maxStep=maxStep, data=initOut$data,
                  verbose = verbose)
out <- fitEMOut$out
Est.points <- fitEMOut$Est.points
phi_fletcher <- out$phi_fletcher
my.design.matrix <- mgcv::model.matrix.gam(out$GamObj)
lengthUniqueDataID <- length(unique(initOut$data$ID))
phi_reml <- out$GamObj$reml.scale
estimateBZOut <- estimateBZ(Posit=initOut$data$Posit,
                            my.design.matrix=my.design.matrix,
                            ncolsZ=ncol(Z), n.k=n.k, verbose = verbose)
Beta.out <- estimateBeta(BZ=estimateBZOut$BZ, BZ.beta=estimateBZOut$BZ.beta,
                         n.k=n.k, Z=Z, out=out, verbose = verbose)
cum_s <- cumsum(n.k)
w_ij <- out$pi.ij * (1 - out$pi.ij)/phi_fletcher
initOut_data <- initOut$data
inputs <- list(w_ij=w_ij, new.par=out$par, new.lambda=out$lambda,
               X=initOut_data$X, Y=initOut_data$Y, my.design.matrix=my.design.matrix,
               gam_smoothMat=out$GamObj$smooth, Z=Z, pred.pi=out$pi.ij,
               p0=p0, p1=p1, disp_est=phi_fletcher, RanEff=RanEff,
               N=lengthUniqueDataID, verbose = verbose)
load(paste0(getwd(), "/data/ref_input_hessianComp.rda"))

H <- hessianComp(w_ij = inputs$w_ij, new.par=inputs$new.par,
                 new.lambda=inputs$new.lambda,  X=inputs$X,
                 Y=inputs$Y, my.design.matrix=inputs$my.design.matrix,
                 gam_smoothMat=inputs$gam_smoothMat, Z=inputs$Z,
                 pred.pi=inputs$pred.pi, p0=inputs$p0,
                 p1=inputs$p1,
                 disp_est=inputs$disp_est, RanEff=inputs$RanEff,
                 N= inputs$N, verbose = FALSE)

test_that("hessianComp calculates the correct Hessian matrix", {
  expect_equal(nrow(H), length(inputs$new.par))
  expect_equal(ncol(H), length(inputs$new.par))
  expect_true(isTRUE(all.equal(solve(-H), BEM.obj$cov1)))
  #expect_true(all(solve(-H) - ref$cov1 <0.000001))
})
