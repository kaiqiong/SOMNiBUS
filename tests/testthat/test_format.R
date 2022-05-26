context("testing parsing functions")

# set seed for reproducibility
set.seed(1234)

M <- matrix(1:9, 3,3)
Cov <- M + 2
colnames(M) <- c("A1", "A2", "A3")
BS1 <- bsseq::BSseq(pos = 1:3, chr = c("chr1", "chr2", "chr1"), M = M, Cov = Cov)
dat1 <- formatFromBSseq(bsseq_dat = BS1, verbose = FALSE)


M <- matrix(1:12, 3,2)
Cov <- M + 2
colnames(M) <- c("A1", "A2")
BS2 <- bsseq::BSseq(pos = 1:3, chr = rep("chr1",3), M = M, Cov = Cov)
dat2 <- formatFromBSseq(bsseq_dat = BS2, verbose = FALSE)


M <- matrix(sample(x = 0:10, size = 12, replace = T), 6,2)
Cov <- M + 2
M[3,2] <- M[6,1] <- Cov[3,2] <- Cov[6,1] <- 0
colnames(M) <- c("A1", "A2")
BS3 <- bsseq::BSseq(pos = 1:6, chr = rep("chr1",6), M = M, Cov = Cov)
dat3 <- formatFromBSseq(bsseq_dat = BS3, verbose = FALSE)


test_that("The length of results matches the expected length", {
  expect_true(length(dat1) == 2)
  expect_true(length(dat2) == 1)
  expect_true(length(dat3) == 1)
})

# check that 'bad data' are filtered out
expect_true(object = nrow(dat3$chr1) == 10)

infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov.gz",
                      package = "bsseq")
dat4 <- formatFromBismark(infile, verbose = FALSE)

# check the format
expect_true(sum(!colnames(dat4[[1]]) %in% c("ID","Position","Meth_Counts","Total_Counts")) == 0)

indir <- system.file("extdata/ev_bt2_tab", package = "SOMNiBUS")
dat5 <- formatFromBSmooth(indir, verbose = FALSE)

# check the format
expect_true(sum(!colnames(dat5[[1]]) %in% c("ID","Position","Meth_Counts","Total_Counts")) == 0)

# change that the additional arguments are passed to the bsseq function
dat5 <- formatFromBSmooth(indir, sampleNames = "test", verbose = FALSE)
expect_true(unique(dat5[[1]]$ID) == "test")

