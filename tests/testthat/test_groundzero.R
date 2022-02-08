context("testing that we recover the same values over and over with a given dataset")

data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# BEM.obj <- binomRegMethModel(
#   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#   Quasi = FALSE, RanEff = FALSE, verbose = FALSE)
load(paste0(getwd(), "/data/BEM.obj.rda"))
# unnecessary test 
# default parameters to be used
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

out_1 <- binomRegMethModel(
  data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon,
  epsilon.lambda = epsilon.lambda, maxStep = maxStep, 
  Quasi = Quasi, RanEff = RanEff, verbose = FALSE
)

test_that("all values in reference output, saved on disk, equals actual values coming from binomRegMethModel", {
  expect_true(isTRUE(all.equal(out_1, BEM.obj, tolerance = 1e-2)))
})
