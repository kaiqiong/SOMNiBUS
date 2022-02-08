context("testing splitData")
require("IRanges")
# set seed for reproducibility
set.seed(1234)

# create input
blocks <- c(1:100,201:250,1251:1350,1361:1380,1401:1500)

my_data <- data.frame(Meth_Counts = sample.int(n = 20, size = length(blocks), 
                                               replace = TRUE), 
                      Total_Counts = sample.int(n = 70, size = length(blocks), 
                                                replace = TRUE),
                      Position = blocks, ID = 1)

#----#

# small gap, no minimal CpGs
small_gap_no_min_cpgs <- splitDataByRegion(dat = my_data, gap = 8, min.cpgs = 1, verbose = FALSE)

# medium gap, no minimal CpGs
medium_gap_no_min_cpgs <- splitDataByRegion(dat = my_data, gap = 80, min.cpgs = 1, verbose = FALSE)

# big gap, no minimal CpGs
big_gap_no_min_cpgs <- splitDataByRegion(dat = my_data, gap = 500, min.cpgs = 1, verbose = FALSE)

# no gap, no minimal CpGs
no_gap_no_min_cpgs <- splitDataByRegion(dat = my_data, gap = 1e+6, min.cpgs = 1, verbose = FALSE)

test_that("The length of results matches the expected length", {
  expect_true(length(small_gap_no_min_cpgs) == 5)
  expect_true(length(medium_gap_no_min_cpgs) == 3)
  expect_true(length(big_gap_no_min_cpgs) == 2)
  expect_true(length(no_gap_no_min_cpgs) == 1)
})

#----#

# small gap, small minimal CpGs
small_gap_small_min_cpgs <- splitDataByRegion(dat = my_data, gap = 8, 
                                              min.cpgs = 50, verbose = FALSE)

# medium gap, small minimal CpGs
medium_gap_small_min_cpgs <- splitDataByRegion(dat = my_data, gap = 80, 
                                               min.cpgs = 50, verbose = FALSE)

# big gap, small minimal CpGs
big_gap_small_min_cpgs <- splitDataByRegion(dat = my_data, gap = 500,
                                            min.cpgs = 50, verbose = FALSE)

# no gap, small minimal CpGs
no_gap_small_min_cpgs <- splitDataByRegion(dat = my_data, gap = 1e+6,
                                           min.cpgs = 50, verbose = FALSE)

test_that("The length of results matches the expected length", {
  expect_true(length(small_gap_small_min_cpgs) == 4)
  expect_true(length(medium_gap_small_min_cpgs) == 3)
  expect_true(length(big_gap_small_min_cpgs) == 2)
  expect_true(length(no_gap_small_min_cpgs) == 1)
})

#----#

# small gap, medium minimal CpGs
small_gap_medium_min_cpgs <- splitDataByRegion(dat = my_data, gap = 8, 
                                               min.cpgs = 100, verbose = FALSE)

# medium gap, medium minimal CpGs
medium_gap_medium_min_cpgs <- splitDataByRegion(dat = my_data, gap = 80, 
                                                min.cpgs = 100, verbose = FALSE)

# big gap, medium minimal CpGs
big_gap_medium_min_cpgs <- splitDataByRegion(dat = my_data, gap = 500, 
                                             min.cpgs = 100, verbose = FALSE)

# no gap, medium minimal CpGs
no_gap_medium_min_cpgs <- splitDataByRegion(dat = my_data, gap = 1e+6,
                                            min.cpgs = 100, verbose = FALSE)

test_that("The length of results matches the expected length", {
  expect_true(length(small_gap_medium_min_cpgs) == 3)
  expect_true(length(medium_gap_medium_min_cpgs) == 2)
  expect_true(length(big_gap_medium_min_cpgs) == 2)
  expect_true(length(no_gap_medium_min_cpgs) == 1)
})

#----#

# small gap, big minimal CpGs
small_gap_big_min_cpgs <- splitDataByRegion(dat = my_data, gap = 8, 
                                            min.cpgs = 200, verbose = FALSE)

# medium gap, big minimal CpGs
medium_gap_big_min_cpgs <- splitDataByRegion(dat = my_data, gap = 50, 
                                             min.cpgs = 200, verbose = FALSE)

# big gap, big minimal CpGs
big_gap_big_min_cpgs <- splitDataByRegion(dat = my_data, gap = 500,
                                          min.cpgs = 200, verbose = FALSE)

# no gap, big minimal CpGs
no_gap_big_min_cpgs <- splitDataByRegion(dat = my_data, gap = 1e+6,
                                         min.cpgs = 200, verbose = FALSE)

test_that("The length of results matches the expected length", {
  expect_true(length(small_gap_big_min_cpgs) == 0)
  expect_true(length(medium_gap_big_min_cpgs) == 1)
  expect_true(length(big_gap_big_min_cpgs) == 1)
  expect_true(length(no_gap_big_min_cpgs) == 1)
})

#----#

# small gap, small minimal CpGs, small maximal CpGs
small_gap_small_min_max_cpgs <- splitDataByRegion(dat = my_data, gap = 8, 
                                                  min.cpgs = 20, max.cpgs = 50, verbose = FALSE)

# medium gap, small minimal CpGs, small maximal CpGs
medium_gap_small_min_max_cpgs <- splitDataByRegion(dat = my_data, gap = 80, 
                                                   min.cpgs = 20, max.cpgs = 50, verbose = FALSE)

# big gap, small minimal CpGs, small maximal CpGs
big_gap_small_min_max_cpgs <- splitDataByRegion(dat = my_data, gap = 500,
                                                min.cpgs = 20, max.cpgs = 50, verbose = FALSE)

# no gap, small minimal CpGs, small maximal CpGs
no_gap_small_min_max_cpgs <- splitDataByRegion(dat = my_data, gap = 1e+6,
                                               min.cpgs = 20, max.cpgs = 50, verbose = FALSE)

test_that("The length of results matches the expected length", {
  expect_true(length(small_gap_small_min_max_cpgs) == 2)
  expect_true(length(medium_gap_small_min_max_cpgs) == 1)
  expect_true(length(big_gap_small_min_max_cpgs) == 0)
  expect_true(length(no_gap_small_min_max_cpgs) == 0)
})

#----#


# small gap, small minimal CpGs, medium maximal CpGs
small_gap_small_min_medium_max_cpgs <- splitDataByRegion(dat = my_data, gap = 8, 
                                                         min.cpgs = 20, max.cpgs = 100, verbose = FALSE)

# medium gap, small minimal CpGs, medium maximal CpGs
medium_gap_small_min_medium_max_cpgs <- splitDataByRegion(dat = my_data, gap = 80, 
                                                          min.cpgs = 20, max.cpgs = 100, verbose = FALSE)

# big gap, small minimal CpGs, medium maximal CpGs
big_gap_small_min_medium_max_cpgs <- splitDataByRegion(dat = my_data, gap = 500,
                                                       min.cpgs = 20, max.cpgs = 100, verbose = FALSE)

# no gap, small minimal CpGs, medium maximal CpGs
no_gap_small_min_medium_max_cpgs <- splitDataByRegion(dat = my_data, gap = 1e+6,
                                                      min.cpgs = 20, max.cpgs = 100, verbose = FALSE)

test_that("The length of results matches the expected length", {
  expect_true(length(small_gap_small_min_medium_max_cpgs) == 5)
  expect_true(length(medium_gap_small_min_medium_max_cpgs) == 2)
  expect_true(length(big_gap_small_min_medium_max_cpgs) == 0)
  expect_true(length(no_gap_small_min_medium_max_cpgs) == 0)
})

#----#

# small gap, small minimal CpGs, big maximal CpGs
small_gap_small_min_big_max_cpgs <- splitDataByRegion(dat = my_data, gap = 8, 
                                                      min.cpgs = 20, max.cpgs = 300, verbose = FALSE)

# medium gap, small minimal CpGs, big maximal CpGs
medium_gap_small_min_big_max_cpgs <- splitDataByRegion(dat = my_data, gap = 80, 
                                                       min.cpgs = 20, max.cpgs = 300, verbose = FALSE)

# big gap, small minimal CpGs, big maximal CpGs
big_gap_small_min_big_max_cpgs <- splitDataByRegion(dat = my_data, gap = 500,
                                                    min.cpgs = 20, max.cpgs = 300, verbose = FALSE)

# no gap, small minimal CpGs, big maximal CpGs
no_gap_small_min_big_max_cpgs <- splitDataByRegion(dat = my_data, gap = 1e+6,
                                                   min.cpgs = 20, max.cpgs = 300, verbose = FALSE)

test_that("The length of results matches the expected length", {
  expect_true(length(small_gap_small_min_big_max_cpgs) == 5)
  expect_true(length(medium_gap_small_min_big_max_cpgs) == 3)
  expect_true(length(big_gap_small_min_big_max_cpgs) == 2)
  expect_true(length(no_gap_small_min_big_max_cpgs) == 0)
})

#----#

#--test region dropping--# 
m1 <- capture_messages(splitDataByRegion(dat = my_data, gap = 80, 
                                         min.cpgs = 100))
m2 <- capture_messages(splitDataByRegion(dat = my_data, gap = 500, 
                                         min.cpgs = 200))
m3 <- capture_messages(splitDataByRegion(dat = my_data, gap = 8,
                                         min.cpgs = 50))
m4 <- capture_messages(splitDataByRegion(dat = my_data, gap = 500,
                                         min.cpgs = 20, max.cpgs = 200))
m5 <- capture_messages(splitDataByRegion(dat = my_data, gap = 500,
                                         min.cpgs = 20, max.cpgs = 200,
                                         verbose = FALSE))

test_that(desc = "Test region dropping",{
  expect_match(object = m1[2], regexp = ".*201-250.*50.*100")
  expect_match(object = m2[2], regexp = ".*1-250.*50.*200")
  expect_match(object = m3[2], regexp = ".*1361-1380.*20.*50")
  expect_match(object = m4[2], regexp = ".*1251-1500.*220.*200")
  expect_true(length(m5) == 0)
})

#--test incorrect values--# 
test_that(desc = "Test incorrect values",{
  expect_warning(object = splitDataByRegion(dat = my_data, gap = -10, 
                                            min.cpgs = 50), 
                 regexp = "\\'gap\\' should be a positive integer\\. We will use the default value 1e\\+6 \\(1Mb\\).")
  expect_warning(object = splitDataByRegion(dat = my_data, gap = "10", 
                                            min.cpgs = 50), 
                 regexp = "\\'gap\\' should be a positive integer\\. We will use the default value 1e\\+6 \\(1Mb\\).")
  expect_warning(object = splitDataByRegion(dat = my_data, gap = 10,
                                            min.cpgs = -50), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByRegion(dat = my_data, gap = 10, 
                                            min.cpgs = "50"), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByRegion(dat = my_data, gap = 10,
                                            max.cpgs = -50), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByRegion(dat = my_data, gap = 10, 
                                            max.cpgs = "50"), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByRegion(dat = my_data, gap = 10, 
                                            max.cpgs = 50, min.cpgs = 2000), 
                 regexp = "\\'max.cpgs\\' should be a higher than \\'min.cpgs\\'\\. We will use the default values: min.cpgs=50 and max.cpgs=2000.")
})
