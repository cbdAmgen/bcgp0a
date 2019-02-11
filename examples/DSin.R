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

xTrain <- matrix(seq(0, 3, length.out = 15), ncol = 1)
yTrain <- DSin(xTrain) + rnorm(length(xTrain), 0, 0.02)
xPred <- matrix(sort(c(xTrain, seq(min(xTrain), max(xTrain), 0.005))), ncol = 1)

noise <- FALSE
priors <- createPriors(xTrain, noise = noise)
inits <- createInits(xTrain, priors, chains = 1)
chains <- 1

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 3, numAdapt = 1000,
            burnin = 1000, nmcmc = 5000, chains <- 1, cores = 1,
            noise = noise)

lapply(fit@sim$samples[[1]], mean)
plot(xTrain, colMeans(fit@sim$samples[[1]]$V), type = 'l')
