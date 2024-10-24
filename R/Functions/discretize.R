# Discretize  Distributions
# Epidemia-style discretization of gamma
epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, y)
  pmf[1] <- pgamma(1.5, alpha, rate = beta)
  for (i in 2:y) {
    pmf[i] <- pgamma(i + .5, alpha, rate = beta) -
      pgamma(i - .5, alpha, rate = beta)
  }
  return(pmf)
}

zero_epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, (y + 1))
  pmf[1] <- pgamma(0.5, alpha, rate = beta)
  for (i in 2:(y + 1)) {
    pmf[i] <- pgamma(i - 1 + .5, alpha, rate = beta) -
      pgamma(i - 1 - .5, alpha, rate = beta)
  }
  return(pmf)
}

epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, y)
  pmf[1] <- sdprisk::phypoexp(1.5, rates)
  for (i in 2:y) {
    pmf[i] <- sdprisk::phypoexp(i + .5, rates) -
      sdprisk::phypoexp(i - .5, rates)
  }
  return(pmf)
}

epidemia_lognormal <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- plnorm(1.5, meanlog = params[1], sdlog = params[2])
  for (i in 2:y) {
    pmf[i] <- plnorm(i + .5, meanlog = params[1], sdlog = params[2]) -
      plnorm(i - .5, meanlog = params[1], sdlog = params[2])
  }
  return(pmf)
}

epidemia_weibull <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- pweibull(1.5, shape = params[1], scale = params[2])
  for (i in 2:y) {
    pmf[i] <- pweibull(i + .5, shape = params[1], scale = params[2]) -
      pweibull(i - .5, shape = params[1], scale = params[2])
  }
  return(pmf)
}

zero_epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, (y + 1))
  pmf[1] <- sdprisk::phypoexp(0.5, rates)
  for (i in 2:(y + 1)) {
    pmf[i] <- sdprisk::phypoexp(i - 1 + .5, rates) -
      sdprisk::phypoexp(i - 1 - .5, rates)
  }
  return(pmf)
}
