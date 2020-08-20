context("testing hessianComp")

path_ref_results <- paste(paste(getwd(), "/data/", sep = ""), "ref_results.RDS",
                                                sep = ""
)
ref <- readRDS(path_ref_results)

path_ref_input_hessianComp <- paste(paste(getwd(), "/data/", sep = ""), "ref_input_hessianComp.RDS",
                                                sep = ""
)
ref_input_hessianComp <- readRDS(path_ref_input_hessianComp)

H <- hessianComp(w_ij=w_ij, new.par=out$par,
                 new.lambda=out$lambda,  X=fitGamOut$data$X, Y=fitGamOut$data$Y, my.design.matrix=my.design.matrix, gam.int=out$GamObj, Z=Z, pred.pi=out$pi.ij, p0=p0, p1=p1,
                 disp_est=phi_fletcher, RanEff=RanEff, N=lengthUniqueDataID)

test_that("hessianComp calculates the correct Hessian matrix", {
    expect_equal(nrow(H), length(out$par))
    expect_equal(ncol(H), length(out$par))
    expect_true(isTRUE(all.equal(solve(-H), ref$cov1)))
    #expect_true(all(solve(-H) - ref$cov1 <0.000001))
})

