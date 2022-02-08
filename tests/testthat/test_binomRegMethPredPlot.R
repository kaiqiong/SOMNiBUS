context("testing plot_binomRegMethPred")

load(paste0(getwd(), "/data/BEM.obj.rda"))

# create newdata
pos <- BEM.obj$uni.pos
my.p <- length(pos)
newdata <- expand.grid(pos, c(0, 1), c(0, 1))
colnames(newdata) <- c("Position", "T_cell", "RA")
my.pred <- binomRegMethModelPred(BEM.obj, newdata, type = "link.scale")

# prepare data for binomRegMethPredPlot
newdata$group <- ""
newdata[(newdata$RA == 0 & newdata$T_cell == 0),]$group <- "CTRL MONO"
newdata[(newdata$RA == 0 & newdata$T_cell == 1),]$group <- "CTRL TCELL"
newdata[(newdata$RA == 1 & newdata$T_cell == 0),]$group <- "RA MONO"
newdata[(newdata$RA == 1 & newdata$T_cell == 1),]$group <- "RA TCELL"
pred <- cbind(newdata, Pred = my.pred)

p.plot <- binomRegMethPredPlot(pred = pred, pred.type = "link.scale", 
                               pred.col = "Pred", verbose = FALSE)
test_that("binomRegMethPredPlot return ggplot object", {
  expect_is(object = p.plot, class = "ggplot")
})

#--test incorrect values--# 
test_that(desc = "Test incorrect values",{
  expect_warning(object = binomRegMethPredPlot(pred = pred, 
                                               pred.type = "logit", 
                                               pred.col = "Pred", 
                                               verbose = FALSE), 
                 regexp = "Unknown prediction type. We will used the type \\'proportion\\'.")
  expect_warning(object = binomRegMethPredPlot(pred = pred, 
                                               pred.col = "Pred", 
                                               group.col = "group",
                                               style = list(),
                                               verbose = FALSE), 
                 regexp = "Ignoring the custom style because of discrepancies between \\'style\\' list and the group defined in the \\'pred\\' data.frame.")
  expect_error(object = binomRegMethPredPlot(pred = pred, 
                                             pred.type = "link.scale", 
                                             pred.col = "predictions", 
                                             verbose = FALSE), 
               regexp = "The 'pred.col' column name is missing from the 'pred' data.frame")
})

#--test that the plot is saved--#
test_that(desc = "Test incorrect values",{
  filename <- tempfile(fileext = ".png")
  # check that the file doesn't already exist
  expect_false(object = file.exists(filename))
  p.plot <- binomRegMethPredPlot(pred = pred, pred.type = "link.scale", 
                                 pred.col = "Pred", 
                                 save = filename,
                                 verbose = FALSE)
  # check the file has been created by binomRegMethPredPlot
  expect_true(object = file.exists(filename))
  
  # test the handling od unknown format
  filename <- tempfile()
  # check that the file doesn't already exist
  expect_false(object = file.exists(filename))
  expect_warning(object = binomRegMethPredPlot(pred = pred, 
                                               pred.type = "link.scale", 
                                               pred.col = "Pred", 
                                               save = filename, 
                                               verbose = FALSE), 
                 regexp = paste0("Unknown graphics format for \\'",filename,
                                 "\\'. The plot will be saved as a .pdf by default"))
  # check the file has been created by binomRegMethPredPlot
  expect_true(object = file.exists(paste0(filename,".pdf")))
})
