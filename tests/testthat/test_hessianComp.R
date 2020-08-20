context("testing hessianComp")

path_ref_results <- paste(paste(getwd(), "/data/", sep = ""), "ref_results.RDS",
                                                sep = ""
)
ref <- readRDS(path_ref_results)

path_ref_input_hessianComp <- paste(paste(getwd(), "/data/", sep = ""), "ref_input_hessianComp.RDS",
                                                sep = ""
)
ref_input_hessianComp <- readRDS(path_ref_input_hessianComp)

H <- hessianComp(w_ij = ref_input_hessianComp$w_ij, new.par=ref_input_hessianComp$new.par,
                 new.lambda=ref_input_hessianComp$new.lambda,  X=ref_input_hessianComp$X,
                 Y=ref_input_hessianComp$Y, my.design.matrix=ref_input_hessianComp$my.design.matrix,
                 gam.int=ref_input_hessianComp$gam.int, Z=ref_input_hessianComp$Z,
                 pred.pi=ref_input_hessianComp$pred.pi, p0=ref_input_hessianComp$p0,
                 p1=ref_input_hessianComp$p1,
                 disp_est=ref_input_hessianComp$disp_est, RanEff=ref_input_hessianComp$RanEff,
                 N= ref_input_hessianComp$N)

test_that("hessianComp calculates the correct Hessian matrix", {
    expect_equal(nrow(H), length(ref_input_hessianComp$new.par))
    expect_equal(ncol(H), length(ref_input_hessianComp$new.par))
    expect_true(isTRUE(all.equal(solve(-H), ref$cov1)))
    #expect_true(all(solve(-H) - ref$cov1 <0.000001))
})

