context("testing pred_BSMethEM")
p <- 0.003307034 # False positive rate (float)
tp <- 0.9 # True positive rate (float)

data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
out_BSMethEM <- BSMethEM(data = RAdat.f, n.k = rep(5, 3), p0 = tp, p1 = fp, epsilon = 10^(-6), epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE)

out_pred_BSMethEM <- pred_BSMethEM(out_BSMethEM)
test_that("the output retrieved from pred_BSMethEM have the same number of elements than the number of samples times the number of genomic positions", {
  expect_equal(length(out_pred_BSMethEM), length(RAdat.f$Total_Counts))
})
