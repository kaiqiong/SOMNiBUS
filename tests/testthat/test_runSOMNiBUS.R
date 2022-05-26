context("testing runSOMNiBUS")
source("helper.R", local = TRUE)

# load example data.frames
data("RAdat")
data("RAdat2")

set.seed(1234)
# reduce the size of the dataset
# ids <- sample(x = 1:length(unique(RAdat$ID)), size = 20, replace = FALSE)
ids <- c(28,16,22,37,9,5,42,4,34,41,26,6,15,14,20,30,24,36,33,21)
RAdat <- RAdat[RAdat$ID %in% ids,]
RAdat2 <- RAdat2[RAdat2$ID %in% ids,]

# pre-filtered data
RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])

# malformed data.frames
RAdat1 <- RAdat3 <- RAdat
RAdat1$Meth_Counts <- NULL
RAdat3$T_cell <- RAdat3$RA <- NULL

# missing mandatory column
expect_error(object = runSOMNiBUS(dat=RAdat1, 
                                  min.cpgs = 5,
                                  n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE),
             regexp = 'Please make sure object \"dat\" have columns named as \"Meth_Counts\", \"Total_Counts\", \"Position\" and \"ID\"'
)

# missing mandatory covariates
expect_error(object = runSOMNiBUS(dat=RAdat3, 
                                  min.cpgs = 5,
                                  n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE),
             regexp = 'The input data should contain at least one covariate'
)

# info: bad data filtered
expect_message(object = runSOMNiBUS(dat=RAdat, 
                                    min.cpgs = 5,
                                    n.k = rep(5,3), p0 = 0.003, p1 = 0.9),
               regexp = "Remove .* rows with Total_Counts equal to 0."
)


# warning: bad split method name, warning for the selection of default approach
expect_warning(object = runSOMNiBUS(dat=RAdat.f, 
                                    split = list(approach = "gap"),
                                    min.cpgs = 5,
                                    n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE), 
               regexp = 'Unknown partitioning method. We will used the default approach "region".'
)

# test the effect of the internal filter
testthat::expect_equal(object = runSOMNiBUS(dat=RAdat, 
                                            split = list(approach = "region", gap = 1e6), 
                                            min.cpgs = 5,
                                            n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE), 
                       expected = runSOMNiBUS(dat=RAdat.f, 
                                              split = list(approach = "region", gap = 1e6), 
                                              min.cpgs = 5,
                                              n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE)
)

bed_file <- paste(paste(getwd(), "/data/", sep = ""),
                  "test_splitDataByBED.bed", sep = "")

test_that(desc = "Test incorrect values while calling partitioning functions",{
  skip_if_run_bitbucket_pipeline()
  expect_warning(object = runSOMNiBUS(dat = RAdat.f,
                                      split = list(approach = "region", gap = -10),
                                      min.cpgs = 50, n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE),
                 regexp = "\\'gap\\' should be a positive integer\\. We will use the default value 1e\\+6 \\(1Mb\\).")
  expect_warning(object = runSOMNiBUS(dat = RAdat.f,
                                      split = list(approach = "density",
                                                   gap = 100,
                                                   window.size = 100),
                                      min.cpgs = 50, n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE),
                 regexp = "\\'gap\\' should be lesser than the 'window.size' value\\. We will use half of window size \\(50\\).")
  expect_error(object = runSOMNiBUS(dat = RAdat.f, split = list(approach = "bed",
                                                                chr = "CHR",
                                                                gap = 0),
                                    min.cpgs = 50, n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE),
               regexp = "\\'chr\\' should be a character vector which length equals the number of rows in \\'dat\\'.")
  expect_error(object = runSOMNiBUS(dat = RAdat.f, split = list(approach = "gene",
                                                                chr = rep("chr4", nrow(RAdat.f)),
                                                                gap = -1,
                                                                types =  "exon",
                                                                organism = "mouse",
                                                                build = "mm10"),
                                    min.cpgs = 50, n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE),
               regexp = "Please install the following package\\(s\\) in order to perform the requested partiniong: \\'TxDb.Mmusculus.UCSC.mm10.knownGene\\', \\'org.Mm.eg.db\\'.")
  expect_error(object = runSOMNiBUS(dat = RAdat.f, split = list(approach = "chromatin",
                                                                chr = rep("chr4", nrow(RAdat.f)),
                                                                gap = -1,
                                                                cell.line = "hmec",
                                                                states = "Promoter"),
                                    min.cpgs = 50, n.k = rep(5,3), p0 = 0.003, p1 = 0.9, verbose = FALSE),
               regexp = "None of the requested states are supported\\. Choose among these supported states: \\'all\\', \\'ActivePromoter\\', \\'Heterochrom\\', \\'Insulator\\', \\'PoisedPromoter\\', \\'RepetitiveCNV\\', \\'Repressed\\', \\'StrongEnhancer\\', \\'TxnElongation\\', \\'TxnTransition\\', \\'WeakEnhancer\\', \\'WeakPromoter\\', \\'WeakTxn\\'.")
  
})

test_that(desc = "Test that the wrapper function and the partitioning + binomRegMethModel work the same",{
  skip_if_run_bitbucket_pipeline()
  
  # create dataset
  RAdat4 <- RAdat2[, c("Meth_Counts","Total_Counts","Position","ID","ACPA4")]
  
  # partitioning by region => 2 regions
  region_gap_100 <- splitDataByRegion(dat = RAdat4, gap = 100, min.cpgs = 50, verbose = FALSE)
  
  n.k_dim1 <- max(1, round(length(unique(region_gap_100[[1]]$Position))/20))
  region_gap_100_r1 <- binomRegMethModel(data = region_gap_100[[1]],
                                         n.k = rep(n.k_dim1,2), p0 = 0.003,
                                         p1 = 0.9, verbose = FALSE)
  
  n.k_dim1 <- max(1, round(length(unique(region_gap_100[[2]]$Position))/20))
  region_gap_100_r2 <- binomRegMethModel(data = region_gap_100[[2]],
                                         n.k = rep(n.k_dim1,2), p0 = 0.003,
                                         p1 = 0.9, verbose = FALSE)
  
  wrapper_region_gap_100 <- runSOMNiBUS(dat = RAdat4,
                                        split = list(approach = "region", gap = 100),
                                        min.cpgs = 50, n.k = rep(5,3), p0 = 0.003,
                                        p1 = 0.9, verbose = FALSE)
  
  testthat::expect_equal(object = region_gap_100_r1,
                         expected = wrapper_region_gap_100[[1]]
  )
  
  testthat::expect_equal(object = region_gap_100_r2,
                         expected = wrapper_region_gap_100[[2]]
  )
  
  # partitioning by density => 1 region
  density_default <- splitDataByDensity(dat = RAdat4, verbose = FALSE)
  
  n.k_dim1 <- max(1, round(length(unique(density_default[[1]]$Position))/20))
  density_r1 <- binomRegMethModel(data = density_default[[1]],
                                  n.k = rep(n.k_dim1,2), p0 = 0.003,
                                  p1 = 0.9, verbose = FALSE)
  
  wrapper_density <- runSOMNiBUS(dat = RAdat4,
                                 split = list(approach = "density"),
                                 n.k = rep(6,2), p0 = 0.003,
                                 p1 = 0.9, verbose = FALSE)
  
  testthat::expect_equal(object = density_r1,
                         expected = wrapper_density[[1]]
  )
  
  # partitioning by chromatin => 1 region
  # a single region so n.k shouldn't be modified
  chromatin_default <- splitDataByChromatin(dat = RAdat4,
                                            chr = rep("chr4", nrow(RAdat4)), cell.line = "hmec", states = "all")
  chromatin_r1 <- binomRegMethModel(data = chromatin_default[[1]],
                                    n.k = rep(10,2), p0 = 0.003,
                                    p1 = 0.9, verbose = FALSE)
  
  wrapper_chromatin <- runSOMNiBUS(dat = RAdat4,
                                   split = list(approach = "chromatin", chr = rep("chr4", nrow(RAdat4)), cell.line = "hmec", states = "all"),
                                   n.k = rep(10,2), p0 = 0.003,
                                   p1 = 0.9, verbose = FALSE)
  
  testthat::expect_equal(object = chromatin_r1,
                         expected = wrapper_chromatin[[1]]
  )
})