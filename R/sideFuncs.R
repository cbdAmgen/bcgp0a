initFunc <- function(initList, priors, x){

  d <- ncol(x)
  n <- nrow(x)
  initReturn <- vector("list")
  initReturn$beta0 <- rnorm(1, 0, 1)
  initReturn$w <- priors$w$lower + rbeta(1, priors$w$alpha, priors$w$beta)*
    (priors$w$upper - priors$w$lower)
  initReturn$rhoG <- rbeta(d, priors$rhoG$alpha, priors$rhoG$beta)
  initReturn$rhoL <- initReturn$rhoG * rbeta(d, priors$rhoL$alpha, priors$rhoL$beta)
  initReturn$sig2eps <- max(.Machine$double.eps,
                            rgamma(1, shape = priors$sig2eps$alpha,
                                   scale = priors$sig2eps$beta))
  initReturn$muV <- rnorm(1, priors$muV$betaV, sqrt(priors$muV$sig2))
  initReturn$rhoV <- rbeta(d, priors$rhoV$alpha, priors$rhoV$beta)
  initReturn$sig2V <- 1/rgamma(1, priors$sig2V$alpha, scale = priors$sig2V$beta)
  K <- initReturn$sig2V * getCorMat(x,initReturn$rhoV) + 1e-10*diag(n)
  initReturn$V <- exp(MASS::mvrnorm(1, initReturn$muV*rep(1, n), K))
  return(initReturn)

}
