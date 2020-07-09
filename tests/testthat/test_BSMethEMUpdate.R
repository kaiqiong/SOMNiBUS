context("testing BSMethEMupdate")

path_ref_input_BSMethEMUpdate <- paste(paste(getwd(), "/data/", sep = ""), "ref_input_BSMethEMUpdate.RDS", sep = "")
ref_input_BSMethEMUpdate <- readRDS(path_ref_input_BSMethEMUpdate)

out_BSMethEMUpdate <- BSMethEMUpdate(data = ref_input_BSMethEMUpdate$data, pi.ij = ref_input_BSMethEMUpdate$old.pi.ij, p0 = ref_input_BSMethEMUpdate$p0, p1 = ref_input_BSMethEMUpdate$p1, n.k = ref_input_BSMethEMUpdate$n.k, binom.link = ref_input_BSMethEMUpdate$binom.link, method = ref_input_BSMethEMUpdate$method, Z = ref_input_BSMethEMUpdate$Z, my.covar.fm = ref_input_BSMethEMUpdate$my.covar.fm, Quasi = ref_input_BSMethEMUpdate$Quasi, scale = ref_input_BSMethEMUpdate$scale)

test_that("the output corresponds to BSMethEMupdate method description", {
  expect_false(length(out_BSMethEMUpdate) == 4)
})
