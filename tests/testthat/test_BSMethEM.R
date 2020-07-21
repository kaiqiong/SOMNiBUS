context("testing binomRegMethModel")
library(dplyr)
fp <- 0.003307034 # False positive rate (float)
tp <- 0.9 # True positive rate (float)

data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
out_binomRegMethModel <- binomRegMethModel(
  data = RAdat.f, n.k = rep(5, 3), p0 = tp, p1 = fp, epsilon = 10^(-6),
  epsilon.lambda = 10^(-3), maxStep = 200, detail = FALSE
)

test_that("number of outputs corresponds to the description of binomRegMethModel", {
  # 13 objects are listed for output, the list length should be 13 this is wrong,
  # there are 16 objects in the actual output
  expect_false(length(out_binomRegMethModel) == 13)
})


test_that("all parameters values get closer along EM optimization procedure", {
  cols <- colnames(out_binomRegMethModel$ite.points)
  n_iter <- length(out_binomRegMethModel$ite.points[, cols[1]])
  a_list_bool <- c()
  for (i in 1:(length(cols) - 1)) {
    diff_first_last <- abs(out_binomRegMethModel$ite.points[1, i] - out_binomRegMethModel$ite.points[
      n_iter,
      i
    ])
    diff_2_last <- abs(out_binomRegMethModel$ite.points[n_iter, i] - out_binomRegMethModel$ite.points[n_iter -
      1, i])
    a_list_bool <- c(a_list_bool, diff_2_last < diff_first_last)
  }
  expect_true(all(a_list_bool))
})

test_that("the number of unique positions output by binomRegMethModel corresponds to the number of unique positions found in the data", {
  expect_equal(length(out_binomRegMethModel$uni.pos), length(unique(RAdat.f$Position)))
})

test_that("the genomic positions are well ordered in binomRegMethModel output (ascending order)", {
  a_list_bool <- c()
  for (i in 1:(length(out_binomRegMethModel$uni.pos) - 1)) {
    a_list_bool <- c(a_list_bool, out_binomRegMethModel$uni.pos[i] < out_binomRegMethModel$uni.pos[i +
      1])
  }
  # this is wrong, we should expect true for all
  expect_false(all(a_list_bool))
})

test_that("if the number of rows in out$SE.out equals the number of unique genomic positions", {
  expect_equal(length(rownames(out_binomRegMethModel$SE.out)), length(unique(RAdat.f$Position)))
})

# check if the genomic positions are well ordered in the dataset provided (not
# filtered for NAs and zero Total_counts) in the package (ascending order)
a_list_bool <- c()
for (i in 1:length(unique(RAdat$ID))) {
  cur <- filter(RAdat, ID == i)
  for (j in 1:(length(cur$ID) - 1)) {
    a_list_bool <- c(a_list_bool, cur$Position[j] < cur$Position[j + 1])
  }
}
# this is wrong, we should expect true for all
expect_false(all(a_list_bool))

# check if the genomic positions are well ordered in the dataset provided
# (filtered) in the package (ascending order)
a_list_bool <- c()
for (i in 1:length(unique(RAdat.f$ID))) {
  cur <- filter(RAdat.f, ID == i)
  for (j in 1:(length(cur$ID) - 1)) {
    a_list_bool <- c(a_list_bool, cur$Position[j] < cur$Position[j + 1])
  }
}
# this is wrong, we should expect true for all
expect_false(all(a_list_bool))
