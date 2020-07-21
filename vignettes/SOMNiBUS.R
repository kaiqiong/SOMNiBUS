## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  ROOT_PACKAGE_PATH <- paste(getwd(), "/", sep = "")
#  devtools::document(ROOT_PACKAGE_PATH)
#  devtools::load_all(ROOT_PACKAGE_PATH)

## ----dependencies, warning=FALSE, message=FALSE-------------------------------
library(SOMNiBUS)

## -----------------------------------------------------------------------------
data(RAdat)

## -----------------------------------------------------------------------------
head(RAdat)

## -----------------------------------------------------------------------------
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])

## -----------------------------------------------------------------------------
out <- binomRegMethModel(data = RAdat.f, n.k = rep(5, 3), p0 = 0.003, p1 = 0.9, Quasi = FALSE, RanEff = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  out.ctype <- binomRegMethModel(data = RAdat.f, n.k = rep(5, 3), p0 = 0.003, p1 = 0.9, covs = "T_cell")

## -----------------------------------------------------------------------------
out$reg.out

## ---- fig.cap="The estimates (solid red lines) and 95% pointwise confidence intervals (dashed red lines) of the intercept, the smooth effect of cell type and RA on methylation levels. ", fig.height= 4, fig.width=9----
binomRegMethModelPlot(out)

## ---- fig.cap="The estimates (solid red lines) and 95% pointwise confidence intervals (dashed red lines) of the intercept, the smooth effect of cell type and RA on methylation levels. (Same ranges of Y axis.)", fig.height= 4, fig.width=9----
binomRegMethModelPlot(out, same.range = TRUE)

## -----------------------------------------------------------------------------
pos <- out$uni.pos
my.p <- length(pos)
newdata <- expand.grid(pos, c(0, 1), c(0, 1))
colnames(newdata) <- c("Position", "T_cell", "RA")

## -----------------------------------------------------------------------------
my.pred <- binomRegMethModelPred(out, newdata, type = "link.scale")

## ----  fig.height= 6, fig.width=6, fig.cap="The predicted methylation levels in the logit scale for the 4 groups of samples with different disease and cell type status."----
plot(pos[order(pos)], (my.pred[(newdata$RA == 0 & newdata$T_cell == 0)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "Predicted methylation levels (in logit scale)", col = "blue", main = "Logit scale", ylim = c(min(my.pred), max(my.pred)), lwd = 2
)
lines(pos[order(pos)], (my.pred[(newdata$RA == 0 & newdata$T_cell == 1)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "predicted", col = "green", lwd = 2
)
lines(pos[order(pos)], (my.pred[(newdata$RA == 1 & newdata$T_cell == 0)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "predicted", col = "red", lwd = 2
)
lines(pos[order(pos)], (my.pred[(newdata$RA == 1 & newdata$T_cell == 1)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "predicted", col = "black", lwd = 2
)
legend("top", c("RA MONO", "RA TCELL", "CTRL MONO", "CTRL TCELL"),
  fill = c("red", "black", "blue", "green"),
  title = "Disease and Cell Type", bty = "n", cex = 0.8
)

## ---- fig.height= 6, fig.width=6 , fig.cap="The predicted methylation levels  proportion scale (right) for the 4 groups of samples with different disease and cell type status."----
my.pred <- binomRegMethModelPred(out, newdata, type = "proportion")
plot(pos[order(pos)], (my.pred[(newdata$RA == 0 & newdata$T_cell == 0)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "Predicted methylation levels (in logit scale)", col = "blue", main = "Proportion scale", ylim = c(min(my.pred), max(my.pred)), lwd = 2
)
lines(pos[order(pos)], (my.pred[(newdata$RA == 0 & newdata$T_cell == 1)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "predicted", col = "green", lwd = 2
)
lines(pos[order(pos)], (my.pred[(newdata$RA == 1 & newdata$T_cell == 0)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "predicted", col = "red", lwd = 2
)
lines(pos[order(pos)], (my.pred[(newdata$RA == 1 & newdata$T_cell == 1)])[order(pos)],
  type = "l", xlab = "Position",
  ylab = "predicted", col = "black", lwd = 2
)
legend("top", c("RA MONO", "RA TCELL", "CTRL MONO", "CTRL TCELL"),
  fill = c("red", "black", "blue", "green"),
  title = "Disease and Cell Type", bty = "n", cex = 0.8
)

## -----------------------------------------------------------------------------
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])

## -----------------------------------------------------------------------------
time0 <- Sys.time()
out_quasi <- binomRegMethModel(data = RAdat.f, n.k = rep(5, 3), p0 = 0.003, p1 = 0.9, Quasi = TRUE, RanEff = FALSE)
print(Sys.time() - time0)
out_quasi$reg.out

## -----------------------------------------------------------------------------
data("RAdat2")
head(RAdat2)

## -----------------------------------------------------------------------------
nCpG <- length(unique(RAdat2$Position))
time0 <- Sys.time()
out_RA2 <- binomRegMethModel(data = RAdat2, n.k = c(5, rep(max(3, round(nCpG / 20)), 9)), p0 = 0.003, p1 = 0.9, Quasi = TRUE, RanEff = FALSE)
print(Sys.time() - time0)

## -----------------------------------------------------------------------------
out_RA2$reg.out

## -----------------------------------------------------------------------------
binomRegMethModelPlot(out_RA2, mfrow = c(2, 5))

## -----------------------------------------------------------------------------
time0 <- Sys.time()
out_RA2_noError <- binomRegMethModel(data = RAdat2, n.k = c(5, rep(max(3, round(nCpG / 20)), 9)), p0 = 0, p1 = 1, Quasi = TRUE, RanEff = FALSE)
print(Sys.time() - time0)

## -----------------------------------------------------------------------------
out_RA2_noError$reg.out

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()
