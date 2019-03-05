context("vector times matrix inverse times vector")

set.seed(112358)
n <- 10
x <- matrix(seq(0, 1, length.out = n), ncol = 1)
A <- getCorMat(x, rho = 0.6)
y <- MASS::mvrnorm(1, rep(0, n), A)
onesN <- rep(1, n)
x1 <- list(onesN, onesN)
x2 <- list(onesN, y)

oneAinvOne <- onesN %*% solve(A) %*% onesN
oneAinvY <- onesN %*% solve(A) %*% y


testthat::test_that("x1 times matrix inverse times x2 is calculated correctly", {

  expect_equal(x1Ainvx2(x1, A, x2), c(oneAinvOne, oneAinvY))

})
