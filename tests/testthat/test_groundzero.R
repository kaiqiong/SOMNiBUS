context("testing that we recover the same values over and over with a given dataset")

#seed_1 <- 1

data(RAdat)
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])

# default parameters to be used
epsilon <- 10^(-6)
epsilon.lambda <- 10^(-3)
p0 <- 0.003307034
p1 <- 0.9
maxStep <- 200
detail <- FALSE
n.k <- rep(5, 3)

path_ref_data <- paste(paste(getwd(), "/data/", sep = ""), "ref_resultsV2.RDS", sep = "")
ref <- readRDS(path_ref_data)

#set.seed(seed_1)
time.0 <- Sys.time()
out_1 <- binomRegMethModel(
  data = RAdat.f, n.k = n.k, p0 = p0, p1 = p1, epsilon = epsilon,
  epsilon.lambda = epsilon.lambda, maxStep = maxStep, detail = detail
)
print(Sys.time()-time.0)
# this was used to save output object when we started working on that directory:
# 'saveRDS(out_1, path_ref_data)'

test_that("all values in reference output, saved on disk, equals actual values coming from binomRegMethModel", {
  expect_true(isTRUE(all.equal(out_1, ref,tolerance = 1e-3)))
})
