initFunc <- function(initList, priors, xTrain){

  d <- ncol(xTrain)
  n <- nrow(xTrain)
  initReturn <- vector("list")
  initReturn$mu <- rnorm(1, 0, 1)
  initReturn$w <- priors$w$a + rbeta(1, priors$w$alpha, priors$w$beta)*(priors$w$b - priors$w$a)
  initReturn$rhoG <- rbeta(d, priors$rhoG$alpha, priors$rhoG$beta)
  initReturn$rhoL <- initReturn$rhoG * rbeta(d, priors$rhoL$alpha, priors$rhoL$beta)
  initReturn$sig2eps <- max(.Machine$double.eps,
                            rgamma(1, shape = priors$sig2eps$alpha,
                                   scale = priors$sig2eps$beta))
  initReturn$muV <- rnorm(1, priors$muV$mu, sqrt(priors$muV$sig2))
  initReturn$rhoV <- rbeta(d, priors$rhoV$alpha, priors$rhoV$beta)
  initReturn$sig2K <- 1/rgamma(1, priors$sig2K$alpha, scale = priors$sig2K$beta)
  K <- initReturn$sig2K * getCorMat(xTrain,initReturn$rhoV) + 1e-10*diag(n)
  initReturn$V <- exp(MASS::mvrnorm(1, initReturn$muV*rep(1, n), K))
  return(initReturn)

}
