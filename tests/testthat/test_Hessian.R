context("testing Hessian")

test_that("we don't know much about the parameters passed and what is computed here", {
  # this is wrong, we should define a proper test here
  expect_error(Hessian(1))
})
