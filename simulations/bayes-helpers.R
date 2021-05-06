# ==============================================================================
# Acknowledgement: 
# All functions obtained from the BFDA package (https://github.com/nicebread/BFDA)
# With the exception of the EV function (to calculate the mean of the posterior distribution),
# Which was generously given to me by Angelika Stefan. Thank you, Angelika!

# ==============================================================================
# Functions obtained from the BFDA package
# ==============================================================================
# These are the functions for t-tests with informed priors
# ==============================================================================
# see also https://arxiv.org/abs/1704.02479 for the formulae

#' @import hypergeo

# helper functions for the computation of the Bayes factor with informed priors

A <- function(t, n, nu, mu.delta, g) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = mu.delta^2*t^2/
                             (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
}

B <- function(t, n, nu, mu.delta, g) {
  
  out <- mu.delta*t/sqrt(1/2*(1/n + g)*((1 + n*g)*nu + t^2)) *
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2)) *
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = mu.delta^2*t^2/
                               (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
  return(out)
  
}


C <- function(delta, t, n, nu) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = n*t^2*delta^2/(2*(nu + t^2))))
  
}

D <- function(delta, t, n, nu) {
  
  out <- t*delta*sqrt(2*n/(nu + t^2))*
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2))*
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = n*t^2*delta^2/(2*(nu + t^2))))
  
  return(out)
  
}

term_normalprior <- function(t, n, nu, mu.delta, g) {
  
  (1 + n*g)^(-1/2) * exp(-mu.delta^2/(2*(1/n + g))) *
    (1 + t^2/(nu*(1 + n*g)))^(-(nu + 1)/2) *
    (A(t, n, nu, mu.delta, g) + B(t, n, nu, mu.delta, g))
  
}

posterior_normal_tmp <- function(delta, t, n1, n2 = NULL,
                                 independentSamples = FALSE, prior.mean,
                                 prior.variance,
                                 rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, n1*n2/(n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.mean
  g <- prior.variance
  
  numerator <- exp(-neff/2*delta^2)*(1 + t^2/nu)^(-(nu + 1)/2)*
    (C(delta, t, neff, nu) + D(delta, t, neff, nu))*
    dnorm(delta, mu.delta, sqrt(g))
  
  denominator <- term_normalprior(t = t, n = neff, nu = nu,
                                  mu.delta = mu.delta, g = g)
  
  out <- numerator/denominator
  
  if ( is.na(out))
    out <- 0
  
  return(out)
  
}

posterior_normal <- Vectorize(posterior_normal_tmp, "delta")

# ==============================================================================
# Function obtained from Angelika Stefan
# ==============================================================================
EV = function(x, t, n1, n2, independentSamples = TRUE, prior.mean = 0.5, prior.variance = 0.3) {
  
  x * posterior_normal(x, t, n1, n2, independentSamples, prior.mean, prior.variance)
  
}