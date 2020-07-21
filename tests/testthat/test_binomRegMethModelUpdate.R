context("testing binomRegMethModelupdate")

path_ref_input_binomRegMethModelUpdate <- paste(paste(getwd(), "/data/", sep = ""), "ref_input_BSMethEMUpdate.RDS",
  sep = ""
)
ref_input_binomRegMethModelUpdate <- readRDS(path_ref_input_binomRegMethModelUpdate)

out_binomRegMethModelUpdate <- binomRegMethModelUpdate(
  data = ref_input_binomRegMethModelUpdate$data, pi.ij = ref_input_binomRegMethModelUpdate$old.pi.ij,
  p0 = ref_input_binomRegMethModelUpdate$p0, p1 = ref_input_binomRegMethModelUpdate$p1, n.k = ref_input_binomRegMethModelUpdate$n.k,
  binom.link = ref_input_binomRegMethModelUpdate$binom.link, method = ref_input_binomRegMethModelUpdate$method,
  Z = ref_input_binomRegMethModelUpdate$Z, my.covar.fm = ref_input_binomRegMethModelUpdate$my.covar.fm,
  Quasi = ref_input_binomRegMethModelUpdate$Quasi, scale = ref_input_binomRegMethModelUpdate$scale
)

test_that("the output corresponds to binomRegMethModelupdate method description", {
  expect_false(length(out_binomRegMethModelUpdate) == 4)
})
