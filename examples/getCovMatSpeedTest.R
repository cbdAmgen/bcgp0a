getCovMat1 <- function(V, R, sig2eps){

  C <- diag(V)^0.5 %*% R %*% diag(V)^0.5 + sig2eps*diag(1, length(V))
  return(C)

}

getCovMat2 <- function(V, R, sig2eps){

  V <- diag(V)
  C <- V^0.5 %*% R %*% V^0.5 + sig2eps*diag(1, nrow(V))
  return(C)

}

getCovMat3 <- function(V, R, sig2eps){

  rootV <- diag(sqrt(V))
  C <- rootV %*% R %*% rootV + sig2eps*diag(1, nrow(rootV))
  return(C)

}

getCovMat4 <- function(V, R, sig2eps){

  rootV <- Matrix::Diagonal(x = sqrt(V))
  C <- rootV %*% R %*% rootV + sig2eps*diag(1, nrow(rootV))
  return(C)

}

getCovMat5 <- function(V, R, sig2eps){

  rootV <- sqrt(V)
  C <- t(t(rootV*R) * rootV) + sig2eps*diag(1, length(V))
  return(C)

}

getCovMat6 <- function(V, R, sig2eps){

  rootV <- sqrt(V)
  C <- t(rootV * t(rootV*R)) + sig2eps*diag(1, length(V))
  return(C)

}

getCovMat7 <- function(V, R, sig2eps){

  rootV <- sqrt(V)
  C <- t(t(rootV*R) * rootV) + diag(sig2eps, length(V))
  return(C)

}


n <- 5000
d <- 2
x <- matrix(runif(n * d), nrow = n, ncol = d)
rho <- runif(d, 0, 1)
R <- getCorMat(x, rho)
sig2eps <- 0.01

V <- exp(rnorm(n))

# my_check <- function(values) {
#   all(sapply(values[-1], function(x) identical(values[[1]], x)))
# }
#
# microbenchmark::microbenchmark(getCovMat1(V, R, sig2eps),
#                                getCovMat2(V, R, sig2eps),
#                                getCovMat3(V, R, sig2eps),
#                                getCovMat4(V, R, sig2eps),
#                                getCovMat5(V, R, sig2eps),
#                                getCovMat6(V, R, sig2eps),
#                                getCovMat7(V, R, sig2eps))

microbenchmark::microbenchmark(getCovMat4(V, R, sig2eps),
                               getCovMat7(V, R, sig2eps),
                               times = 10)
