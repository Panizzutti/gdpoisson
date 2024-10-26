% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdp_functions.R
\name{logLikgd}
\alias{logLikgd}
\title{Log-Likelihood Function for a sample of the GD-Poisson Distribution}
\usage{
logLikgd(k, lambda, theta)
}
\arguments{
\item{k}{A vector of non-negative integers representing the observed count data.}

\item{lambda}{A positive parameter representing the mean of the distribution.}

\item{theta}{A positive parameter representing the dispersion of the distribution.}
}
\value{
A single numeric value representing the sum of the log-likelihoods for the provided data and parameters.
}
\description{
Calculates the log-likelihood for observed count data under the GD-Poisson distribution. This function is capable of handling both scalar and vector inputs for the observed counts (\code{k}), mean parameters (\code{lambda}), and dispersion parameters (\code{theta}). If any of these parameters are vectors, they must either be of the same length or have a length of 1, in which case the single value is used for all evaluations.
}
\details{
The GD-Poisson (Poisson-Negative Zero) distribution is a generalization of the Poisson distribution that introduces arbitrary underdispersion. The log-likelihood is computed by first calculating the probability mass function (PMF) for each observed count \eqn{k} using the \code{dgdpois_scalar} function. To ensure numerical stability and avoid taking the logarithm of zero or \code{NA} values, any PMF values that are less than or equal to zero or \code{NA} are replaced with a very small positive number (\eqn{\epsilon = 1 \times 10^{-64}}).
}
\examples{
# Single set of parameters
k_obs <- 3
lambda <- 2
theta <- 1
log_likelihood <- logLikgd(k_obs, lambda, theta)
print(log_likelihood)

# Vector of observed counts with scalar lambda and theta
k_obs <- 0:5
log_likelihood <- logLikgd(k_obs, lambda = 2, theta = 1)
print(log_likelihood)

# Vectors of the same length for k, lambda, and theta
k_obs <- 0:5
lambda <- c(2, 2, 3, 3, 4, 4)
theta <- c(1, 1.1, 1.2, 1.3, 1.4, 1.5)
log_likelihood <- logLikgd(k_obs, lambda, theta)
print(log_likelihood)

# Example with invalid input lengths (should throw an error)
tryCatch({
  k_obs <- 0:5
  lambda <- c(2, 3)  # Length 2
  theta <- 1          # Length 1
  log_likelihood <- logLikgd(k_obs, lambda, theta)
}, error = function(e) {
  message("Error: ", e$message)
})

}
