context("testing plot_binomRegMethModel")

### creation of the BEM.obj
# data(RAdat)
# RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
# BEM.obj <- binomRegMethModel(
#   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#   Quasi = FALSE, RanEff = FALSE, verbose = FALSE)
load(paste0(getwd(), "/data/BEM.obj.rda"))

p.plot <- binomRegMethModelPlot(BEM.obj, same.range = FALSE, 
                                verbose = FALSE)
test_that("binomRegMethModelPlot return ggplot object", {
  expect_is(object = p.plot, class = "ggplot")
})

#--test incorrect values--# 
test_that(desc = "Test incorrect values",{
  expect_warning(object = binomRegMethModelPlot(BEM.obj, 
                                                same.range = FALSE, 
                                                covs = c("T_cell", "RA","test"),
                                                verbose = FALSE), 
                 regexp = "Some covariates are missing from the BEM.obj object: \\'test\\'")
  
  expect_error(object = binomRegMethModelPlot(BEM.obj, 
                                              same.range = FALSE, 
                                              covs = c("T_cell2", "RA2","test"),
                                              verbose = FALSE), 
               regexp = "None of the requested covariates are included in the BEM.obj object")
})

#--test that the plot is saved--#
test_that(desc = "Test incorrect values",{
  filename <- tempfile(fileext = ".png")
  # check that the file doesn't already exist
  expect_false(object = file.exists(filename))
  p.plot <- binomRegMethModelPlot(BEM.obj, same.range = FALSE, 
                                  save = filename, verbose = FALSE)
  # check the file has been created by binomRegMethModelPlot
  expect_true(object = file.exists(filename))
  
  # test the handling od unknown format
  filename <- tempfile()
  # check that the file doesn't already exist
  expect_false(object = file.exists(filename))
  expect_warning(object = binomRegMethModelPlot(BEM.obj, same.range = FALSE, 
                                                save = filename, verbose = FALSE), 
                 regexp = paste0("Unknown graphics format for \\'",filename,
                                 "\\'. The plot will be saved as a .pdf by default"))
  # check the file has been created by binomRegMethModelPlot
  expect_true(object = file.exists(paste0(filename,".pdf")))
})
