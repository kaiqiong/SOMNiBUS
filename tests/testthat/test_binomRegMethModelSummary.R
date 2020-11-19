context("testing binomRegMethModelSummary")

path_ref_results <- paste(paste(getwd(), "/data/", sep = ""), "ref_results.RDS",
                          sep = ""
)
ref <- readRDS(path_ref_results)

path_ref_input_binomRegMethModelSummary <- paste(paste(getwd(), "/data/", sep = ""), "ref_input_binomRegMethModelSummary.RDS",
                                    sep = ""
)
ref_input_binomRegMethModelSummary <- readRDS(path_ref_input_binomRegMethModelSummary)




s.table <- binomRegMethModelSummary( ref_input_binomRegMethModelSummary$GamObj, ref_input_binomRegMethModelSummary$var.cov.alpha,
                                     ref_input_binomRegMethModelSummary$new.par ,ref_input_binomRegMethModelSummary$edf.out,
                                     ref_input_binomRegMethModelSummary$edf1.out, ref_input_binomRegMethModelSummary$X_d,
                                     ref_input_binomRegMethModelSummary$resi_df, ref_input_binomRegMethModelSummary$Quasi,
                                     ref_input_binomRegMethModelSummary$scale, ref_input_binomRegMethModelSummary$RanEff,
                                     ref_input_binomRegMethModelSummary$Z)


test_that("binomRegMethModelSummary provides the correct regional p-value summary", {
    expect_equal(nrow(s.table), 4)
    expect_equal(ncol(s.table), 3)
    expect_true(isTRUE(all.equal(s.table, ref$reg.out, tolerance = 1e-3)))
})
