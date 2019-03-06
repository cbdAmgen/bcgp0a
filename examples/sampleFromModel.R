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
V <- as.vector(dbeta(xTrain, 6, 2) + 0.15)
# V <- exp(MASS::mvrnorm(1, rep(muV, length(xTrain)), K))
# V <- rep(1, length(xTrain))
C <- getCovMat(V, R, sig2eps)

yTrain <- MASS::mvrnorm(1, rep(beta0, length(xTrain)), C)
plot(xTrain, yTrain, type = 'l')

chains <- 3


priors <- createPriors(xTrain, noise = noise)
inits <- createInits(xTrain, priors, chains = chains)

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 2, numAdapt = 200,
            burnin = 100, nmcmc = 200, chains = chains, cores = 1,
            noise = noise)


lapply(fit@sim$samples[[1]], mean)
# par(mfrow = c(chains, 1))
blah <- matrix(0, ncol = chains, nrow = n)
for(i in 1:chains) blah[,i] <- colMeans(fit@sim$samples[[i]]$V)
plot(xTrain, rep(1, length(xTrain)), type = 'l', ylab = "Mean V", col = "white",
     ylim = c(0, max(blah)))
for(i in 1:chains) lines(xTrain, colMeans(fit@sim$samples[[i]]$V), type = 'l',
                        ylab = "Mean V", col = i)

mu <- matrix(0, ncol = chains, nrow = length(fit@sim$samples[[1]]$beta0))
for(i in 1:chains) mu[,i] <- fit@sim$samples[[i]]$beta0
plot(fit@sim$samples[[1]]$beta0, type = 'l', ylab = "beta0", col = 1,
     ylim = c(min(mu), max(mu)))
for(i in 2:chains) lines(fit@sim$samples[[i]]$beta0, type = 'l',
                         ylab = "beta0", col = i)

# hist(fit@sim$samples[[1]]$lp__)
