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
yTrain <- DSin(xTrain)
xPred <- matrix(sort(c(xTrain, seq(min(xTrain), max(xTrain), 0.005))), ncol = 1)
priors <- createPriors(xTrain, noise = FALSE)
inits <- createInits(xTrain, priors, chains = 1)

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 3, numAdapt = 1000,
            burnin = 1000, nmcmc = 5000, chains = 1, cores = 1,
            noise = FALSE)
