x <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
y <- x[, 1] + sin(x[, 2])
priors <- createPriors(x, noise = FALSE)
blah <- bcgp(x, y, priors = priors)
blah2 <- bcgp(x, y, priors = "default")
blah3 <- bcgp(x, y, priors = "priors")

inits <- createInits(x, priors = priors)
blah <- bcgp(x, y, priors = priors, inits = inits)
blah2 <- bcgp(x, y, priors = "default", inits = "random")
blah3 <- bcgp(x, y, priors = priors, inits = "default")
blah3 <- bcgp(x, y, priors = priors, inits = createInits(x))

n <- 10
d <- 2
x <- matrix(runif(n * d), nrow = n, ncol = d)
rhoG <- runif(d, 0, 1)
rhoL <- runif(d, 0, rhoG)
G <- getCorMat(x, rhoG)
L <- getCorMat(x, rhoL)
w <- runif(1, 0.5, 1)
R <- combineCorMats(w, G, L)


xTrain <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
yTrain <- xTrain[, 1] + sin(xTrain[, 2])
priors <- createPriors(xTrain, noise = FALSE)
inits <- createInits(xTrain, priors, chains = 1)

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 3, numAdapt = 1000,
            burnin = 1000, nmcmc = 5000, chains = 1, cores = 1,
            noise = FALSE)


################################ Check backsolve for substitute of finding inverse  ########################
n <- 1000
d <- 2
x <- matrix(runif(n * d), nrow = n, ncol = d)
rho <- runif(d, 0, 1)
R <- getCorMat(x, rho)
sig2 <- 0.01
V <- rlnorm(n, -0.1, 0.1)
Sigma <- getCovMat(V, R, sig2)

mu <- 1
y <- rnorm(nrow(Sigma), mu, 1)


yMinusMu <- matrix( c(y - mu), ncol = 1)


m1 <- function(Sigma, y, mu){

  yMinusMu <- matrix( c(y - mu), ncol = 1)
  logLike <- -0.5* t(yMinusMu) %*% solve(Sigma) %*% yMinusMu
  return(logLike)

}

m2 <- function(Sigma, y, mu){

  R <- chol(Sigma)
  yMinusMu <- matrix( c(y - mu), ncol = 1)
  logLike <- -0.5 * t(backsolve(R, yMinusMu, transpose = TRUE))  %*% forwardsolve(t(R), yMinusMu )
  return(logLike)

}

all.equal(m1(Sigma, y, mu), m2(Sigma, y, mu))

microbenchmark::microbenchmark(m1(Sigma, y, mu),
                               m2(Sigma, y, mu),
                               times = 100)

# logDet.R in ~/Documents/bcgp/bcgpr/R
# mvrnormRcpp.cpp in ~/rcppPractice/

