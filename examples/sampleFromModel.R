rm(list = ls())
cat("\014")

n <- 27
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
plot(xTrain, yTrain, type = 'l')

chains <- 1


priors <- createPriors(xTrain, noise = noise)
inits <- createInits(xTrain, priors, chains = chains)

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 5, numAdapt = 500,
            burnin = 100, nmcmc = 2000, chains = chains, cores = 1,
            noise = noise)


lapply(fit@sim$samples[[1]], mean)
plot(xTrain, colMeans(fit@sim$samples[[1]]$V), type = 'l')
