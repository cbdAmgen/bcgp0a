rm(list = ls())
cat("\014")

BJX <- function(x){
  if(!is.matrix(x)){
    stop("x should be a matrix.")
  }
  if(dim(x)[2] != 1){
    stop("x must have exactly one column.")
  }
  y <- sin(30*(x - .9)^4)*cos(2*(x - .9)) + (x - .9)/2
  return(y)
}

xTrain <- matrix(c(seq(0, .4 + .4/11, by = .4/11), seq(0.5, 1, by = 0.5/3)), ncol = 1)
yTrain <- BJX(xTrain)
xPred <- matrix(sort(c(xTrain, seq(min(xTrain), max(xTrain), 0.005))), ncol = 1)

chains <- 1


priors <- createPriors(xTrain, noise = FALSE)
inits <- createInits(xTrain, priors, chains = chains)

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 5, numAdapt = 500,
            burnin = 100, nmcmc = 2000, chains = chains, cores = 1,
            noise = FALSE)

