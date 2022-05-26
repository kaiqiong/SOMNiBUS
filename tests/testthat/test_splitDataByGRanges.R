context("testing splitData by annotations in Genomic Ranges")

# set seed for reproducibility
set.seed(1234)

# create input
blocks <- c(1,2,13,14,15,21,22,28,30)

my_data <- data.frame(Meth_Counts = sample.int(n = 20, size = length(blocks), 
                                               replace = TRUE), 
                      Total_Counts = sample.int(n = 70, size = length(blocks), 
                                                replace = TRUE),
                      Position = blocks, ID = 1)
# my_data
# | Meth_Counts| Total_Counts| Position| ID|
# |-----------:|------------:|--------:|--:|
# |          16|           70|        1|  1|
# |           5|           14|        2|  1|
# |          12|           56|       13|  1|
# |          15|           62|       14|  1|
# |           9|            4|       15|  1|
# |           5|            4|       21|  1|
# |           6|           21|       22|  1|
# |          16|           40|       28|  1|
# |           4|           56|       30|  1|


annots <- data.frame(chr = 1, start = c(1,6,14,18), end = c(3,9,19,29))
# annotations
# | chr| start| end|
# |---:|-----:|---:|
# |   1|     1|   3|
# |   1|     6|   9|
# |   1|    14|  19|
# |   1|    18|  29|

annots_gr <- GenomicRanges::makeGRangesFromDataFrame(df = annots,
                                                     seqnames.field = "chr",
                                                     start.field = "start", 
                                                     end.field = "end")
#----#

#--test incorrect values--# 
test_that(desc = "Test incorrect values",{
  expect_error(object = splitDataByGRanges(dat = my_data, chr = "CHR", gap = 0, 
                                           min.cpgs = 1, annots = annots_gr, verbose = FALSE), 
               regexp = "\\'chr\\' should be a character vector which length equals the number of rows in \\'dat\\'.")
  expect_error(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), gap = 0, 
                                           min.cpgs = 1, annots = "", verbose = FALSE), 
               regexp = "\\'annots\\' should inherit the GRanges class.")
  expect_error(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), gap = 0, 
                                           min.cpgs = 1, annots = annots, verbose = FALSE), 
               regexp = "\\'annots\\' should inherit the GRanges class.")
  expect_error(object = splitDataByGRanges(dat = my_data, verbose = FALSE), 
               regexp = "argument \"chr\" is missing, with no default")
  expect_error(object = splitDataByGRanges(dat = my_data, chr = rep("CHR1", nrow(my_data)), 
                                           gap = 0, min.cpgs = 1, annots = annots_gr, verbose = FALSE), 
               regexp = "Incompatible chromosome name. Please see genomeStyles\\(\\) for supported styles.")
  expect_warning(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                             gap = "0", min.cpgs = 1, 
                                             annots = annots_gr, verbose = FALSE), 
                 regexp = "\\'gap\\' should be an integer >= -1\\. We will use the default value -1.")
  expect_warning(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                             gap = -2, min.cpgs = 1, 
                                             annots = annots_gr, verbose = FALSE), 
                 regexp = "\\'gap\\' should be an integer >= -1\\. We will use the default value -1.")
  expect_warning(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                             gap = 0, min.cpgs = "1", 
                                             annots = annots_gr, verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                             gap = 0, min.cpgs = -1, 
                                             annots = annots_gr, verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                             gap = 0, max.cpgs = "1", 
                                             annots = annots_gr, verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                             gap = 0, max.cpgs = -1, 
                                             annots = annots_gr, verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                             gap = 0, max.cpgs = 1, 
                                             min.cpgs = 50, annots = annots_gr, verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a higher than \\'min.cpgs\\'\\. We will use the default values: min.cpgs=50 and max.cpgs=2000.")
})

#-- test chromosome nomenclature
testthat::expect_equal(object = splitDataByGRanges(dat = my_data, 
                                                   chr = rep("chr1", nrow(my_data)), 
                                                   annots = annots_gr, 
                                                   gap = 0, min.cpgs = 1, verbose = FALSE), 
                       expected = splitDataByGRanges(dat = my_data, 
                                                     chr = rep("1", nrow(my_data)),
                                                     annots = annots_gr, 
                                                     gap = 0, min.cpgs = 1, verbose = FALSE))

#-- min.cpgs = 1, gap = -1 => 3 regions ([1,2],[14,15],[21,22,28]) --#
mc_1_no_gap <- splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), annots = annots_gr, 
                                  gap = -1, min.cpgs = 1, verbose = FALSE)

#-- min.cpgs = 1, gap = 0 => 3 regions ([1,2],[13,14,15],[21,22,28,30]) --#
mc_1_gap_0 <- splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), annots = annots_gr,
                                 gap = 0, min.cpgs = 1, verbose = FALSE)

#-- min.cpgs = 1, gap = 1 => 3 regions ([1,2],[13,14,15,21],[21,22,28,30]) --#
mc_1_gap_1 <- splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), annots = annots_gr, 
                                 gap = 1, min.cpgs = 1, verbose = FALSE)

#-- min.cpgs = 1, gap = 3 => 4 regions ([1,2],[2,13],[13,14,15,21,22],
# [14,15,21,22,28,30]) --#
mc_1_gap_3 <- splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), annots = annots_gr, 
                                 gap = 3, min.cpgs = 1, verbose = FALSE)

### --- ###
#-- min.cpgs = 3, gap = -1 => 1 region ([21,22,28])
# + 2 drops ([1,2],[14,15]) --#
mc_3_no_gap <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                                   annots = annots_gr, gap = -1, 
                                                   min.cpgs = 3, verbose = FALSE))

#-- min.cpgs = 3, gap = 0 => 2 regions ([13,14,15],[21,22,28,30])
# + 1 drop ([1,2]) --#
mc_3_gap_0 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)),
                                                  annots = annots_gr, gap = 0, 
                                                  min.cpgs = 3, verbose = FALSE))

#-- min.cpgs = 3, gap = 1 => 2 regions ([13,14,15,21],[21,22,28,30])
# + 1 drop ([1,2]) --#
mc_3_gap_1 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                                  annots = annots_gr, 
                                                  gap = 1, min.cpgs = 3, verbose = FALSE))

#-- min.cpgs = 3, gap = 3 => 2 regions ([13,14,15,21,22],[14,15,21,22,28,30])
# + 2 drops ([1,2],[2,13]) --#
mc_3_gap_3 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                                  annots = annots_gr, gap = 3, 
                                                  min.cpgs = 3, verbose = FALSE))

### --- ###
#-- min.cpgs = 5, gap = 10 => 3 regions ([1,2,13,14,15],[14,15,21,22,28,30],
# [13,14,15,21,22,28,30]) + 1 drop ([1,2,13,14]) --#
mc_5_gap_10 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                                   annots = annots_gr, gap = 10, 
                                                   min.cpgs = 5, verbose = FALSE))

#-------#

#-- min.cpgs = 3, max.cpgs = 3, gap = 0 => 2 regions ([1,2],[13,14,15])
# + 1 drop ([21,22,28,30]) --#
mc_2_max_3_gap_0 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)),
                                                        annots = annots_gr, gap = 0, 
                                                        min.cpgs = 2, max.cpgs = 3, verbose = FALSE))

#-- min.cpgs = 3, max.cpgs = 4, gap = 1 => 2 regions ([13,14,15,21],[21,22,28,30])
# + 3 drop ([1,2]) --#
mc_3_max_4_gap_1 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                                        annots = annots_gr, 
                                                        gap = 1, min.cpgs = 3, 
                                                        max.cpgs = 4, verbose = FALSE))

#-- min.cpgs = 2, max.cpgs = 4, gap = 3 => 2 regions ([1,2],[2,13])
# + 4 drops ([13,14,15,21,22],[14,15,21,22,28,30]) --#
mc_2_max_4_gap_3 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                                        annots = annots_gr, gap = 3, 
                                                        min.cpgs = 2, max.cpgs = 4, verbose = FALSE))

### --- ###
#-- min.cpgs = 5, max.cpgs = 10, gap = 10 => 3 regions ([1,2,13,14,15],[14,15,21,22,28,30],
# [13,14,15,21,22,28,30]) + 1 drop ([1,2,13,14]) --#
mc_5_max_10_gap_10 <- suppressWarnings(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                                          annots = annots_gr, gap = 10, 
                                                          min.cpgs = 5, max.cpgs = 10, verbose = FALSE))

#-------#

test_that("The length of results matches the expected length", {
  expect_true(length(mc_1_no_gap) == 3)
  expect_true(length(mc_1_gap_0) == 3)
  expect_true(length(mc_1_gap_1) == 3)
  expect_true(length(mc_1_gap_3) == 4)
  expect_true(length(mc_3_no_gap) == 1)
  expect_true(length(mc_3_gap_0) == 2)
  expect_true(length(mc_3_gap_1) == 2)
  expect_true(length(mc_3_gap_3) == 2)
  expect_true(length(mc_5_gap_10) == 3)
  expect_true(length(mc_3_gap_0) == 2)
  expect_true(length(mc_3_gap_1) == 2)
  expect_true(length(mc_3_gap_3) == 2)
  expect_true(length(mc_5_gap_10) == 3)
})

#--test region dropping--# 
m1 <- capture_messages(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                          annots = annots_gr, gap = 0, 
                                          min.cpgs = 3))
m2 <- capture_messages(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                          annots = annots_gr, gap = 1, 
                                          min.cpgs = 3))
m3 <- capture_messages(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)), 
                                          annots = annots_gr, gap = 10, 
                                          min.cpgs = 5))
m4 <- capture_messages(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)),
                                          annots = annots_gr, gap = 0, 
                                          min.cpgs = 2, max.cpgs = 3))
m5 <- capture_messages(splitDataByGRanges(dat = my_data, chr = rep("chr1", nrow(my_data)),
                                          annots = annots_gr, gap = 0, 
                                          min.cpgs = 2, max.cpgs = 3,
                                          verbose = FALSE))

test_that(desc = "Test region dropping",{
  expect_match(object = m1[1], regexp = ".*1:1-2.*2.*3")
  expect_match(object = m2[1], regexp = ".*1:1-2.*2.*3")
  expect_match(object = m3[1], regexp = ".*1:1-14.*4.*5")
  expect_match(object = m4[1], regexp = ".*1:21-30.*4.*3")
  expect_true(length(m5) == 0)
})
