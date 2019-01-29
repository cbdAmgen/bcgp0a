## Matrices will be fastest. In particular, extracting from a data.frame is slow. Lists are similar to
## matrices, but roughly 4-5 times slower (DFs are roughly 100 times slower)

matFunc <- function(x){

  # nRow <- nrow(x)
  y <- x[4,]

  toReturn <- y["mu"] * y["w"] + sum(y[startsWith(names(y), "rhoG")]) - sum(y[startsWith(names(y), "rhoL")]) +
    y["sig2eps"] + y["beta0"]*y["sig2V"] + prod(y[startsWith(names(y), "rhoV")])

  return(unname(toReturn))
}

dfFunc <- function(x){

  # nRow <- nrow(x)
  y <- x[4,]

  toReturn <- y["mu"] * y["w"] + sum(y[startsWith(names(y), "rhoG")]) - sum(y[startsWith(names(y), "rhoL")]) +
    y["sig2eps"] + y["beta0"]*y["sig2V"] + prod(y[startsWith(names(y), "rhoV")])

  return(as.numeric(unname(toReturn)))
}


fun1 <- function(lst, n){
  sapply(lst, `[`, n)
}

listFunc <- function(x){

  # nRow <- nrow(x)
  rhoGs <- fun1(x[startsWith(names(x), "rhoG")], 4)
  rhoLs <- fun1(x[startsWith(names(x), "rhoL")], 4)
  rhoVs <- fun1(x[startsWith(names(x), "rhoV")], 4)
  toReturn <- x$mu[4] * x$w[4] + sum(rhoGs) - sum(rhoLs) +
    x$sig2eps[4] + x$beta0[4]*x$sig2V[4] + prod(rhoVs)

  return(toReturn)

}



d <- 3
n <- 10
dMat <- 5 + 3*d

rhoNames <- rep("rho", d)
rhoGNames <- paste0(rhoNames, paste0("G", 1:d))
rhoLNames <- paste0(rhoNames, paste0("L", 1:d))
rhoVNames <- paste0(rhoNames, paste0("V", 1:d))

myMat <- matrix(runif(n*dMat), nrow = n, ncol = dMat)
colnames(myMat) <- c("mu", "w", rhoGNames, rhoLNames, "sig2eps", "beta0", "sig2V", rhoVNames)

myDF <- data.frame(myMat)

myList <- as.list(myDF)

all.equal(matFunc(myMat), dfFunc(myDF), listFunc(myList))



microbenchmark::microbenchmark(m =  matFunc(myMat),
                               d = dfFunc(myDF),
                               l = listFunc(myList))
