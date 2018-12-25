x <- matrix(runif(20, 0, 10), nrow = 10, ncol = 2)
y <- x[, 1] + sin(x[, 2])
priors <- createPriors(x, noise = FALSE)
blah <- bcgp(x, y, priors = priors)
blah2 <- bcgp(x, y, priors = "default")
blah3 <- bcgp(x, y, priors = "priors")

inits <- createInits(x, priors = priors)
blah <- bcgp(x, y, priors = priors, inits = inits)
blah2 <- bcgp(x, y, priors = "default", inits = "random")
blah3 <- bcgp(x, y, priors = priors, inits = "default")
blah3 <- bcgp(x, y, priors = priors, inits = createInits(x))

