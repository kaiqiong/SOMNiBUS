context("testing splitData by density")
require("IRanges")
# set seed for reproducibility
set.seed(1234)

# create input
blocks <- c(1,2,13,14,15,21,22,28,30)

my_data <- data.frame(Meth_Counts = sample.int(n = 20, size = length(blocks), 
                                               replace = TRUE), 
                      Total_Counts = sample.int(n = 70, size = length(blocks), 
                                                replace = TRUE),
                      Position = blocks, ID = 1)

# r <- unique(IRanges::IRanges(start = my_data$Position, 
# end = my_data$Position))
# > knitr::kable(reduce(r))
# | start| end| width|
# |-----:|---:|-----:|
# |     1|   2|     2|
# |    13|  15|     3|
# |    21|  22|     2|
# |    28|  28|     1|
# |    30|  30|     1|

# > gaps(r)
# | start| end| width|
# |-----:|---:|-----:|
# |     3|  12|    10|
# |    16|  20|     5|
# |    23|  27|     5|
# |    29|  29|     1|

#----#

#--test incorrect values--# 

test_that(desc = "Test incorrect values",{
  expect_warning(object = splitDataByDensity(dat = my_data, window.size = "1", 
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'window.size\\' should be a positive integer greater than 10\\. We will use the default value 100.")
  expect_warning(object = splitDataByDensity(dat = my_data, window.size = 1,
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'window.size\\' should be a positive integer greater than 10\\. We will use the default value 100.")
  expect_warning(object = splitDataByDensity(dat = my_data, by = "1",
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'by\\' should be a positive integer\\. We will use the default value 1.")
  expect_warning(object = splitDataByDensity(dat = my_data, by = 0,
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'by\\' should be a positive integer\\. We will use the default value 1.")
  expect_warning(object = splitDataByDensity(dat = my_data, window.size = 10,
                                             gap = 100, min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'gap\\' should be lesser than the \\'window.size\\' value\\. We will use half of window size \\(5\\).")
  expect_warning(object = splitDataByDensity(dat = my_data, min.density = "1",
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'min.density\\' should be a positive integer\\. We will use the default value 5.")
  expect_warning(object = splitDataByDensity(dat = my_data, min.density = 0,
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'min.density\\' should be a positive integer\\. We will use the default value 5.")
  expect_warning(object = splitDataByDensity(dat = my_data, gap = "1",
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'gap\\' should be a positive integer\\. We will use the default value 10.")
  expect_warning(object = splitDataByDensity(dat = my_data, gap = 0,
                                             min.cpgs = 1, verbose = FALSE), 
                 regexp = "\\'gap\\' should be a positive integer\\. We will use the default value 10.")
  expect_warning(object = splitDataByDensity(dat = my_data, gap = 10,
                                            min.cpgs = -50, verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByDensity(dat = my_data, gap = 10, 
                                            min.cpgs = "50", verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByDensity(dat = my_data, gap = 10,
                                            max.cpgs = -50, verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByDensity(dat = my_data, gap = 10, 
                                            max.cpgs = "50", verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByDensity(dat = my_data, gap = 10, 
                                            max.cpgs = 50, min.cpgs = 2000, verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a higher than \\'min.cpgs\\'\\. We will use the default values: min.cpgs=50 and max.cpgs=2000.")
})


#-- window.size = 10, by = 1, min.density = 1,
# min.cpgs = 5, gap = 1 => 1 region --#
ws_10_by_1_md_1_mc_5_gap_1 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 1, 
                                                 min.cpgs = 5, gap = 1, verbose = FALSE)

#-------#

#-- window.size = 10, by = 1, min.density = 2, 
# min.cpgs = 2, gap = 1 => 2 regions --#
ws_10_by_1_md_2_mc_2_gap_1 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 2, 
                                                 min.cpgs = 2, gap = 1, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 2,
# min.cpgs = 5, gap = 1 => 1 region + 1 drop --#
ws_10_by_1_md_2_mc_5_gap_1 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 2, 
                                                 min.cpgs = 5, gap = 1, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 2,
# min.cpgs = 5, gap = 2 => 1 region + 1 drop --#
ws_10_by_1_md_2_mc_5_gap_2 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 2, 
                                                 min.cpgs = 5, gap = 2, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 2,
# min.cpgs = 5, gap = 3 => 1 region [9 cpgs] --#
ws_10_by_1_md_2_mc_5_gap_3 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 2, 
                                                 min.cpgs = 5, gap = 3, verbose = FALSE)

#-------#


#-- window.size = 10, by = 1, min.density = 3,
# min.cpgs = 2, gap = 1 => 2 dropped regions --#
ws_10_by_1_md_3_mc_5_gap_1 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 3, 
                                                 min.cpgs = 5, gap = 1, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 3,
# min.cpgs = 2, gap = 1 => 2 regions --#
ws_10_by_1_md_3_mc_2_gap_1 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 3, 
                                                 min.cpgs = 2, gap = 1, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 3,
# min.cpgs = 2, gap = 2 => 2 regions --#
ws_10_by_1_md_3_mc_5_gap_2 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 3, 
                                                 min.cpgs = 2, gap = 2, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 3,
# min.cpgs = 5, gap = 3 => 1 regions --#
ws_10_by_1_md_3_mc_5_gap_3 <- splitDataByDensity(dat = my_data, 
                                                 window.size = 10, 
                                                 by = 1, min.density = 3, 
                                                 min.cpgs = 5, gap = 3, verbose = FALSE)

#--- testing max.cpgs ---#

#-- window.size = 10, by = 1, min.density = 1,
# min.cpgs = 5, max.cpgs = 10, gap = 1 => 1 region --#
ws_10_by_1_md_1_max_10_gap_1 <- splitDataByDensity(dat = my_data, 
                                                   window.size = 10, 
                                                   by = 1, min.density = 1, 
                                                   min.cpgs = 5, max.cpgs = 10,
                                                   gap = 1, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 2, 
# min.cpgs = 2, max.cpgs = 3, gap = 1 => 1 region (1 drop) --#
ws_10_by_1_md_2_max_3_gap_1 <- splitDataByDensity(dat = my_data, 
                                                  window.size = 10, 
                                                  by = 1, min.density = 2, 
                                                  min.cpgs = 2, max.cpgs = 3,
                                                  gap = 1, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 2,
# min.cpgs = 3, max.cgps = 5. gap = 1 => 0 region (2 drops) --#
ws_10_by_1_md_2_min_3_max_5_gap_1 <- splitDataByDensity(dat = my_data, 
                                                        window.size = 10, 
                                                        by = 1, min.density = 2, 
                                                        min.cpgs = 3, max.cpgs = 5,
                                                        gap = 1, verbose = FALSE)

#-- window.size = 10, by = 1, min.density = 2,
# min.cpgs = 5, max.cpgs = 10, gap = 3 => 1 region [9 cpgs] --#
ws_10_by_1_md_2_min_5_max_10_gap_3 <- splitDataByDensity(dat = my_data, 
                                                         window.size = 10, 
                                                         by = 1, min.density = 2, 
                                                         min.cpgs = 5, max.cpgs = 10,
                                                         gap = 3, verbose = FALSE)

#-------#

test_that("The length of results matches the expected length", {
  expect_true(length(ws_10_by_1_md_1_mc_5_gap_1) == 1)
  expect_true(length(ws_10_by_1_md_2_mc_2_gap_1) == 2)
  expect_true(length(ws_10_by_1_md_2_mc_5_gap_1) == 1)
  expect_true(length(ws_10_by_1_md_2_mc_5_gap_2) == 1)
  expect_true(length(ws_10_by_1_md_2_mc_5_gap_3) == 1)
  expect_true(length(ws_10_by_1_md_3_mc_5_gap_1) == 0)
  expect_true(length(ws_10_by_1_md_3_mc_2_gap_1) == 2)
  expect_true(length(ws_10_by_1_md_3_mc_5_gap_2) == 2)
  expect_true(length(ws_10_by_1_md_3_mc_5_gap_3) == 1)
  expect_true(length(ws_10_by_1_md_1_max_10_gap_1) == 1)
  expect_true(length(ws_10_by_1_md_2_max_3_gap_1) == 1)
  expect_true(length(ws_10_by_1_md_2_min_3_max_5_gap_1) == 0)
  expect_true(length(ws_10_by_1_md_2_min_5_max_10_gap_3) == 1)
})








