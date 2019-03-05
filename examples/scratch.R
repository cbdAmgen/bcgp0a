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
n <- 20
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


m1 <- function(Sigma, yMinusMu){

  # yMinusMu <- matrix( c(y - mu), ncol = 1)
  logLike <- -0.5* t(yMinusMu) %*% solve(Sigma) %*% yMinusMu
  return(logLike)

}

m2 <- function(Sigma, yMinusMu){

  R <- chol(Sigma)
  # yMinusMu <- matrix( c(y - mu), ncol = 1)
  logLike <- -0.5 * t(backsolve(R, yMinusMu, transpose = TRUE))  %*% forwardsolve(t(R), yMinusMu )
  return(logLike)

}

m3 <- function(Sigma, yMinusMu){

  R <- chol(Sigma)
  # yMinusMu <- matrix( c(y - mu), ncol = 1)
  tmp <- forwardsolve(t(R), yMinusMu )
  logLike <- -0.5 * sum(tmp^2)
  return(logLike)

}
all.equal(m1(Sigma, yMinusMu), m2(Sigma, yMinusMu), m3(Sigma, yMinusMu))

microbenchmark::microbenchmark(m1(Sigma, y, mu),
                               m2(Sigma, y, mu),
                               m3(Sigma, y, mu),
                               times = 100)

# logDet.R in ~/Documents/bcgp/bcgpr/R
# mvrnormRcpp.cpp in ~/rcppPractice/



##################### chol(C) vs. svd(C) vs qr.solve(C) #########################
rm(list = ls())
cat("\014")
x <- matrix(seq(0, 3, length.out = 50), ncol = 1)
y <- exp(-2*x) * sin(4*pi*x^2)
x <- matrix(runif(20, 0, 1), nrow = 10, ncol = 2)
y <- x[, 1] + sin(x[, 2])
priors <- createPriors(x, noise = FALSE)
inits <- createInits(x, priors, chains = 4)
numUpdates <- 3
numAdapt <- 500
burnin <- 500
nmcmc <- 5000
chains <- 4
cores <- 1

nTrain <- nrow(x)
d <- ncol(x)
iterations <- numUpdates*numAdapt + burnin + nmcmc
epsV <- 1e-10
tau2 <- 0.08

rhoNames <- rep("rho", d)
rhoGNames <- paste0(rhoNames, paste0("G", 1:d))
rhoLNames <- paste0(rhoNames, paste0("L", 1:d))
rhoVNames <- paste0(rhoNames, paste0("V", 1:d))
rm(rhoNames)

i <- 1
j <- 2

allDraws <- matrix(NA, nrow = iterations, ncol = 5 + 3*d + nTrain)
row1 <- unlist(inits[[i]])
colnames(allDraws) <- names(row1)
allDraws[1, ] <- row1
if(d == 1){
  colnames(allDraws)[startsWith(colnames(allDraws), "rho")] <-
    paste0(colnames(allDraws)[startsWith(colnames(allDraws), "rho")],"1")
}

propWidths <- c(0.4, rep(0.07, d), rep(0.03, d), 1e-3, rep(0.25, d))
names(propWidths) <- c("w", rhoGNames, rhoLNames, "sig2eps", rhoVNames)



G <- getCorMat(x, allDraws[j - 1, rhoGNames])
L <- getCorMat(x, allDraws[j - 1, rhoLNames])
R <- combineCorMats(allDraws[j - 1, "w"], G, L)
C <- getCovMat(allDraws[j - 1, startsWith(colnames(allDraws), "V")], R,
               allDraws[j - 1, "sig2eps"])
K <- allDraws[j - 1, "sig2V"] * getCorMat(x, allDraws[j - 1, rhoVNames]) +
  diag(epsV, nTrain)

allDraws[j, ] = allDraws[j - 1, ]

## get a sample for beta0. Gibbs step
# b <- 15
# cholCR <- chol(C[1:b, 1:b])
# tmpC <- forwardsolve(t(cholCR), rep(1, b))
# tmpC2 <- forwardsolve(t(cholCR), y[1:b])
# oneCinvone <- sum(tmpC^2)
# oneCinvY <- sum(tmpC * tmpC2)

cholCR <- chol(C)
tmpC <- forwardsolve(t(cholCR), rep(1, nTrain))
tmpC2 <- forwardsolve(t(cholCR), y)
oneCinvone <- sum(tmpC^2)
oneCinvY <- sum(tmpC * tmpC2)

oneCinvone
oneCinvY


# svdC <- svd(C[1:b, 1:b])
svdC <- svd(C)
Cinv <- svdC$v %*% diag(1/svdC$d) %*% t(svdC$u)
sum(rep(1, nTrain) * (Cinv %*% rep(1, nTrain)))
sum(rep(1, nTrain) * (Cinv %*% y))
# sum(rep(1, b) * (Cinv %*% rep(1, b)))
# sum(rep(1, b) * (Cinv %*% y[1:b]))

sum(rep(1, nTrain) * qr.solve(C, rep(1, nTrain)))
sum(rep(1, nTrain) * qr.solve(C, y))

rm(list = ls())
cat("\014")

n <- 1500
beta0 <- 0
w <- 0.65
rhoG <- 0.5
rhoL <- 0.25
sig2eps <- 0
muV <- -0.1
sig2V <- 0.01
rhoV <- 0.09
xTrain <- matrix(seq(0, 1, length.out = n), ncol = 1)
noise <- FALSE

G <- getCorMat(xTrain, rhoG)
L <- getCorMat(xTrain, rhoL)
R <- combineCorMats(w, G, L)
K <- sig2V*getCorMat(xTrain, rhoV) + diag(1e-10, length(xTrain))
# V <- exp(MASS::mvrnorm(1, rep(muV, length(xTrain)), K))
V <- rep(1, length(xTrain))
C <- getCovMat(V, R, sig2eps)

determinant(C)
det(C)
exp(determinant(C)$modulus)
logDet(C)

yTrain <- MASS::mvrnorm(1, rep(beta0, length(xTrain)), C)


#################################################################################

rm(list = ls())
cat("\014")

n <- 10
beta0 <- 0
w <- 0.65
rhoG <- 0.6
rhoL <- 0.3
sig2eps <- 0
muV <- -0.1
sig2V <- 0.01
rhoV <- 0.99
xTrain <- matrix(seq(0, 1, length.out = n), ncol = 1)
noise <- FALSE

G <- getCorMat(xTrain, rhoG)
L <- getCorMat(xTrain, rhoL)
R <- combineCorMats(w, G, L)
K <- sig2V*getCorMat(xTrain, rhoV) + diag(1e-10, length(xTrain))
# V <- exp(MASS::mvrnorm(1, rep(muV, length(xTrain)), K))
V <- rep(1, length(xTrain))
C <- getCovMat(V, R, sig2eps)
yTrain <- MASS::mvrnorm(1, rep(beta0, length(xTrain)), C)


x1 <- rep(1, n)
x2 <- yTrain
A <- C

x1Ainvx2 <- function(x1, A, x2){

  cholAR <- try(chol(A), silent = TRUE)
  if(is.matrix(cholAR)){

    tmp1 <- forwardsolve(t(cholAR), x1)
    tmp2 <- forwardsolve(t(cholAR), x2)
    xAinvy <- sum(tmp1 * tmp2)

  }else{

    Ainvy <- try(solve(A, x2), silent = TRUE)
    if(is.numeric(Ainvy)){
      xAinvy <- sum(x1 * Ainvy)
    }else{
      svdA <- svd(A)
      xAinvy <- t(x1) %*% svdA$v %*% diag(1/svdA$d) %*%
        t(svdA$u) %*% x2
    }

  }
  return(xAinvy)
}

x1Ainvx2(x1, A, x1); t(x1) %*% solve(A) %*% x1
x1Ainvx2(x1, A, x2); t(x1) %*% solve(A) %*% yTrain


###############################
x1 <- list(rep(1, n), rep(1,n))
x2 <- list(rep(1, n), yTrain)
A <- C

x1Ainvx2Lists <- function(x1, A, x2){

  cholAR <- try(chol(A), silent = TRUE)
  if(is.matrix(cholAR)){

    tmp1 <- lapply(x1, forwardsolve, l = t(cholAR))
    tmp2 <- lapply(x2, forwardsolve, l = t(cholAR))

    xAinvy <- colSums(mapply(`*`, tmp1, tmp2))

  }else{

    Ainvy <- try(lapply(x2, solve, a = A), silent = TRUE)
    if(all(sapply(Ainvy, is.numeric))){
      xAinvy <- colSums(mapply(`*`, x1, Ainvy))
    }else{
      svdA <- svd(A)
      middle <- svdA$v %*% diag(1/svdA$d) %*% t(svdA$u)
      end <- lapply(x2, `%*%`, middle)
      xAinvy <- colSums(mapply(`*`, x1, end))
    }
  }
  return(xAinvy)
}

x1Ainvx2Lists(x1, A, x2)



# ## Sample for beta0. Gibbs step
# cholCR <- try(chol(C), silent = TRUE)
# if(is.matrix(cholCR)){
#   tmpC <- forwardsolve(t(cholCR), rep(1, nTrain))
#   tmpC2 <- forwardsolve(t(cholCR), y)
#   oneCinvone <- sum(tmpC^2)
#   oneCinvY <- sum(tmpC * tmpC2)
# }else{
#   Cinvone <- try(solve(C, onesNTrain), silent = TRUE)
#   if(is.numeric(Cinvone)){
#     oneCinvone <- sum(onesNTrain * Cinvone)
#   }else{
#     svdC <- svd(C)
#     oneCinvone <- t(onesNTrain) %*% svdC$v %*% diag(1/svdC$d) %*%
#       t(svdC$u) %*% onesNTrain
#   }
#   CinvY <- try(solve(C, y), silent = TRUE)
#   if(is.numeric(CinvY)){
#     oneCinvY <- sum(onesNTrain * CinvY)
#   }else{
#     svdC <- svd(C)
#     oneCinvY <- t(onesNTrain) %*% svdC$v %*% diag(1/svdC$d) %*%
#       t(svdC$u) %*% y
#   }
# }


x1Ainvx2CBD <- function(x1, A, x2){

  cholAR <- try(chol(A), silent = TRUE)
  if(is.matrix(cholAR)){


    tmp1 <- lapply(x1, forwardsolve, l = t(cholAR))
    tmp2 <- lapply(x2, forwardsolve, l = t(cholAR))

    xAinvy <- colSums(mapply(`*`, tmp1, tmp2))

  }else{

    Ainvy <- try(lapply(x2, solve, a = A), silent = TRUE)
    if(all(sapply(Ainvy, is.numeric))){
      xAinvy <- colSums(mapply(`*`, x1, Ainvy))
    }else{
      svdA <- svd(A)
      middle <- svdA$v %*% diag(1/svdA$d) %*% t(svdA$u)
      end <- lapply(x2, `%*%`, middle)
      xAinvy <- colSums(mapply(`*`, x1, end))
    }
  }
  return(xAinvy)
}
