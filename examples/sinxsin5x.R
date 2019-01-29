rm(list = ls())
cat("\014")

sinxsin5x <- function(x){
  if(!is.matrix(x)){
    stop("x should be a matrix.")
  }
  if(dim(x)[2] != 2){
    stop("x must have exactly two columns.")
  }
  y <- sin(x[, 1]) + sin(5*x[, 2]);
  return(y)
}


n = 24;
d = 2;

xTrain <- matrix(runif(n*d), ncol = d, nrow = n);
yTrain <- sinxsin5x(xTrain)
xPred <- expand.grid(seq(0, 1, length.out= 10),
                     seq(0, 1, length.out= 10))

chains <- 3

priors <- createPriors(xTrain, noise = FALSE)
inits <- createInits(xTrain, priors, chains = chains)

fit <- bcgp(x = xTrain, y = yTrain, priors = priors,
            inits = inits, numUpdates = 2, numAdapt = 100,
            burnin = 100, nmcmc = 200, chains = chains, cores = 1,
            noise = FALSE)


