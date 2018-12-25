context("Correlation matrix")

x <- seq(.01, 0.99, length = 5)
mat <- as.matrix(expand.grid(x, x, x))
rho <- c(0.3, 0.65, 0.9)

R <- getCorMat(x = mat, rho = rho)

testthat::test_that("Correlation matrix is calculated correctly", {

  expect_equal(getCorMat(x = mat, rho = rho), R)

})
