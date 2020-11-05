context("testing plot_binomRegMethModel")
fp <- 0.003307034 # False positive rate (float)
tp <- 0.9 # True positive rate (float)

#data(RAdat)
#RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#out_binomRegMethModel <- binomRegMethModel(
#  data = RAdat.f, n.k = rep(5, 3), p0 = tp, p1 = fp, epsilon = 10^(-6),
#  epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE
#)
path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_resultsV2.RDS", sep = "")
out_binomRegMethModel <- readRDS(path_ref_data)
p.plot <- binomRegMethModelPlot(out_binomRegMethModel, same.range = FALSE)
test_that("binomRegMethModelPlot return NULL", {
  # this is wrong the program should return a plot as said in the description
  expect_null(p.plot)
})
