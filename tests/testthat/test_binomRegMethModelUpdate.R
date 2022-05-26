context("testing binomRegMethModelupdate")

### creation of the BEM.obj
# data(RAdat)
# RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# BEM.obj <- binomRegMethModel(
#   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#   Quasi = FALSE, RanEff = FALSE, verbose = FALSE)
load(paste0(getwd(), "/data/BEM.obj.rda"))

### creation of the inputs for binomRegMethModelUpdate
### this code is coming from the binomRegMethModel and fitEM functions 
### just before calling the fitEM function
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
# 
# # coming from binomRegMethModel
# initOut <- binomRegMethModelInit(data=data, covs=covs, n.k=n.k, 
#                                  verbose = verbose)
# Z <- initOut$Z
# fitGamOut <- fitGam(data=initOut$data, Quasi=Quasi, binom.link=binom.link,
#                     method=method, RanEff=RanEff, scale=scale, Z=Z, n.k=n.k,
#                     verbose = verbose)
# phi_fletcher <- phiFletcher(data=initOut$data, Quasi=Quasi, 
#                             reml.scale=reml.scale, scale=scale, 
#                             gam.int=fitGamOut$gam.int, verbose=verbose)
# 
# # coming from fitEM
# old.pi.ij <-fitGamOut$gam.int$fitted.values
# mean_part <- mean(((old.pi.ij-p0)*(p1-old.pi.ij))/(old.pi.ij*(1-old.pi.ij)),
#                   na.rm = TRUE)
# phi_fletcher <- (phi_fletcher-1)/mean_part+1
# 
# inputs <- list(data=initOut$data, 
#                old.pi.ij=fitGamOut$gam.int$fitted.values,
#                p0=p0, p1=p1, n.k=n.k, 
#                binom.link=binom.link, method=method,
#                Z=Z, my.covar.fm=fitGamOut$my.covar.fm, 
#                Quasi=Quasi, scale=phi_fletcher, 
#                reml.scale=reml.scale, verbose = verbose)
load(paste0(getwd(), "/data/ref_input_BSMethEMUpdate.rda"))

out_binomRegMethModelUpdate <- binomRegMethModelUpdate(
  data = inputs$data, pi.ij = inputs$old.pi.ij,
  p0 = inputs$p0, p1 = inputs$p1, n.k = inputs$n.k,
  binom.link = inputs$binom.link, method = inputs$method,
  Z = inputs$Z, my.covar.fm = inputs$my.covar.fm,
  Quasi = inputs$Quasi, scale = inputs$scale, verbose = FALSE
)

test_that("the output corresponds to binomRegMethModelupdate method description", {
  expect_true(length(out_binomRegMethModelUpdate) == 10)
})


