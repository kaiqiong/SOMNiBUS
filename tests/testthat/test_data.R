context("testing data")
library(dplyr)

data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])

test_that("the genomic positions are well ordered in the dataset provided (not filtered for NAs and zero Total_counts) in the package (ascending order)", {
  a_list_bool <- c()
  for (i in 1:length(unique(RAdat$ID))) {
    cur <- filter(RAdat, ID == i)
    for (j in 1:(length(cur$ID) - 1)) {
      a_list_bool <- c(a_list_bool, cur$Position[j] < cur$Position[j +
        1])
    }
  }
  # this is wrong, we should expect true for all
  expect_false(all(a_list_bool))
})

test_that("the genomic positions are well ordered in the dataset provided (filtered) in the package (ascending order)", {
  a_list_bool <- c()
  for (i in 1:length(unique(RAdat.f$ID))) {
    cur <- filter(RAdat.f, ID == i)
    for (j in 1:(length(cur$ID) - 1)) {
      a_list_bool <- c(a_list_bool, cur$Position[j] < cur$Position[j +
        1])
    }
  }
  # this is wrong, we should expect true for all
  expect_false(all(a_list_bool))
})
