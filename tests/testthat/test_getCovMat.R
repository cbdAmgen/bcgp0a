context("Covariance matrix")

xx <- seq(.01, 0.99, length = 5)
X <- as.matrix(expand.grid(xx, xx, xx))
rho <- c(0.3, 0.65, 0.9)

n <- nrow(X)
d <- ncol(X)
R <- matrix(0, ncol = n, nrow = n)

for(i in 1:n){
  for(j in 1:n){
    if(j >= i){
      dist <- 4 * (X[i, ] - X[j, ])
      R[i,j] = prod(rho^(dist^2))
    }else{
      R[i,j] <- R[j,i]
    }
  }
}

sig2 <- 0.01
V <- rlnorm(n, -0.1, 0.1)
C <- diag(V)^0.5 %*% R %*% diag(V)^0.5 + diag(sig2, n)


testthat::test_that("Covariance matrix is calculated correctly", {

  expect_equal(getCovMat(V, R, sig2), C)

})
