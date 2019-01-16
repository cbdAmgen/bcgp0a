context("Correlation matrix")

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


testthat::test_that("Correlation matrix is calculated correctly", {

  expect_equal(getCorMat(x = X, rho = rho), R)

})
