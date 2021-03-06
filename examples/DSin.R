rm(list = ls())
cat("\014")

DSin <- function(x){
  if(!is.matrix(x)){
    stop("x should be a matrix.")
  }
  if(dim(x)[2] != 1){
    stop("x must have exactly one column.")
  }
  y <- exp(-2*x) * sin(4*pi*x^2)
  return(y)
}

xTrain <- matrix(seq(0, 3, length.out = 25), ncol = 1)
yTrain <- DSin(xTrain) + rnorm(length(xTrain), 0, 0.00)
xPred <- matrix(sort(c(xTrain, seq(min(xTrain), max(xTrain), 0.005))), ncol = 1)
chains <- 2

noise <- FALSE
priors <- createPriors(xTrain, noise = noise)
inits <- createInits(xTrain, priors, chains = chains)


system.time({
fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 3, numAdapt = 100,
            burnin = 100, nmcmc = 200, chains = chains, cores = 1,
            noise = noise)
})

lapply(fit@sim$samples[[1]], mean)
plot(xTrain, colMeans(fit@sim$samples[[1]]$V), type = 'l')
lapply(fit@sim$acceptances[[1]], mean)
