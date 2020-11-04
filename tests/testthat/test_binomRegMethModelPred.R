context("testing binomRegMethModelPred")
fp <- 0.003307034 # False positive rate (float)
tp <- 0.9 # True positive rate (float)

data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
out_binomRegMethModel <- binomRegMethModel(
  data = RAdat.f, n.k = rep(5, 3), p0 = tp, p1 = fp, epsilon = 10^(-6),
  epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE
)

out_binomRegMethModelPred <- binomRegMethModelPred(out_binomRegMethModel)
test_that("the output retrieved from binomRegMethModelPred have the same number of elements than the number of samples times the number of genomic positions", {
  expect_equal(length(out_binomRegMethModelPred), length(RAdat.f$Total_Counts))
})
