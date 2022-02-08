context("testing binomRegMethModelPred")

data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])

### creation of the BEM.obj
# data(RAdat)
# RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# BEM.obj <- binomRegMethModel(
#   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#   Quasi = FALSE, RanEff = FALSE, verbose = FALSE)
load(paste0(getwd(), "/data/BEM.obj.rda"))

out_binomRegMethModelPred <- binomRegMethModelPred(BEM.obj, verbose = FALSE)
test_that("the output retrieved from binomRegMethModelPred have the same number of elements than the number of samples times the number of genomic positions", {
  expect_equal(length(out_binomRegMethModelPred), length(RAdat.f$Total_Counts))
})
