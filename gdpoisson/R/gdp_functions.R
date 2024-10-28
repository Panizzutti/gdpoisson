#' Cumulative Distribution Function of the GD-Poisson Distribution
#'
#' Computes the cumulative distribution function (CDF) of the GD-Poisson distribution for a vector of integer values. The parameters `lambda` and `theta` can also be vectors of the same length as `q`, or any parameter can be a single value to be applied across all evaluations.
#'
#' @param q A vector of non-negative integers at which to evaluate the CDF.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#'
#' @return A numeric vector of the same length as \code{q}, containing the CDF evaluated at each element of \code{q}.
#' @details The GD-Poisson distribution is a generalization of the Poisson distribution that introduces arbitrary underdispersion.
#' @examples
#' # Evaluate the CDF at q = 0:10 with lambda = 5 and theta = 4
#' pgdpois(0:10, lambda = 5, theta = 6)
#'
#' @export
pgdpois <- function(q, lambda, theta) {
  # Get the lengths of each input
  lengths <- c(length(q), length(lambda), length(theta))

  # Check if all lengths are the same or if any length is 1
  if (!all(lengths == max(lengths) | lengths == 1)) {
    stop("All inputs must either be of the same length or length 1.")
  }

  # Determine the length to use for vectorization
  max_length <- max(lengths)

  # Repeat values if necessary to match the maximum length
  q <- if (length(q) == 1) rep(q, max_length) else q
  lambda <- if (length(lambda) == 1) rep(lambda, max_length) else lambda
  theta <- if (length(theta) == 1) rep(theta, max_length) else theta

  # Apply the scalar function element-wise
  sapply(seq_along(q), function(i) {
    pgdpois_scalar(q[i], lambda[i], theta[i])
  })

}


# Internal function: Scalar CDF of the GD-Poisson distribution
pgdpois_scalar <- function(k, lambda, theta) {
  if (k < 0 ) {
    return(0)
  }
  if ( k != floor(k)){
    k <- floor(k)
  }
  if (lambda <= 0 || theta <= 0) {
    stop("Parameters 'lambda' and 'theta' must be positive.")
  }



  #original
  kp1 <- k + 1
  if (k == 0) {
    return((kp1 - lambda) * (pgamma(lambda, shape = kp1 / theta, scale = theta, lower.tail = FALSE)) +
             kp1 * (dgamma(lambda / theta, shape = 1 + kp1 / theta, scale = 1)))
  }
  return(
    (kp1 - lambda) * (pgamma(lambda, shape = kp1 / theta, scale = theta, lower.tail = FALSE)) +
      kp1 * (dgamma(lambda / theta, shape = 1 + kp1 / theta, scale = 1)) -
      (k - lambda) * (pgamma(lambda, shape = k / theta, scale = theta, lower.tail = FALSE)) -
      k * (dgamma(lambda / theta, shape = 1 + k / theta, scale = 1))
  )


}

# R/dgdpois.R

#' Probability Mass Function of the GD-Poisson Distribution
#'
#' Computes the probability mass function (PMF) of the GD-Poisson distribution for a vector of integer values. The parameters \code{lambda} and \code{theta} can also be vectors of the same length as \code{k}, or any parameter can be a single value to be applied across all evaluations.
#'
#' @param k A vector of non-negative integers at which to evaluate the PMF.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#'
#' @return A numeric vector of the same length as \code{k}, containing the PMF evaluated at each element of \code{k}.
#' @details The GD-Poisson distribution allows for greater flexibility in modeling count data with over-dispersion compared to the standard Poisson distribution.
#'
#' The function is vectorized, meaning that \code{k}, \code{lambda}, and \code{theta} can be vectors. If any of these parameters are single values (length 1), they will be used for all evaluations. All vector inputs must either be of the same length or have a length of 1.
#'
#' @examples
#' # Evaluate the PMF at k = 0:10 with lambda = 2 and theta = 1
#' dgdpois(0:10, lambda = 2, theta = 1)
#'
#' # Vectorized parameters: different values for theta
#' dgdpois(0:5, lambda = 2, theta = c(1, 1.1, 1.2, 1.3, 1.4, 1.5))
#'
#' # Scalar lambda and theta with vector k
#' dgdpois(0:5, lambda = 2, theta = 1)
#'
#' @export
dgdpois <- function(k, lambda, theta) {
  # Get the lengths of each input
  lengths <- c(length(k), length(lambda), length(theta))

  # Check if all lengths are the same or if any length is 1
  if (!all(lengths == max(lengths) | lengths == 1)) {
    stop("All inputs must either be of the same length or length 1.")
  }

  # Determine the length to use for vectorization
  max_length <- max(lengths)

  # Repeat values if necessary to match the maximum length
  k <- if (length(k) == 1) rep(k, max_length) else k
  lambda <- if (length(lambda) == 1) rep(lambda, max_length) else lambda
  theta <- if (length(theta) == 1) rep(theta, max_length) else theta

  # Apply the scalar function element-wise
  sapply(seq_along(k), function(i) {
    dgdpois_scalar(k[i], lambda[i], theta[i])
  })
}


# Internal function: Scalar PMF of the GD-Poisson distribution
dgdpois_scalar <- function(k, lambda, theta) {
  if (k != floor(k)){
    stop("k must be an integer")
  }

  if (k < 0) {
    return(0)
  }

  return(pgdpois_scalar(k, lambda, theta) - pgdpois_scalar(k - 1, lambda, theta))


  }

#' Random Generation from the GD-Poisson Distribution
#'
#' Generates random samples from the GD-Poisson distribution.
#'
#' @param n Integer specifying the number of random samples to generate.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#'
#' @return An integer vector of length \code{n}, containing random samples from the GD-Poisson distribution.
#' @details The generated random sample follow the GD-Poisson distribution.
#' @examples
#' # Generate 10 random samples with lambda = 2 and theta = 1
#' rgdpois(10, lambda = 2, theta = 1)
#' @export
rgdpois <- function(n, lambda, theta) {
  # n: number of samples to generate
  # lambda, theta: parameters of the distribution
  sample <- integer(n)

  for (i in seq_len(n)) {
    u <- runif(1)  # Uniform random number between 0 and 1
    k <- -1        # Start from k = 0 after increment
    repeat {
      k <- k + 1
      if (pgdpois_scalar(k, lambda, theta) >= u) {
        sample[i] <- k
        break
      }
    }
  }
  return(sample)
}

# R/qgdpois.R

#' Quantile Function of the GD-Poisson Distribution
#'
#' Computes the quantile function of the GD-Poisson distribution for a vector of probabilities. The parameters \code{lambda} and \code{theta} can also be vectors of the same length as \code{p}, or any parameter can be a single value to be applied across all evaluations.
#'
#' @param p A vector of probabilities (each between 0 and 1) at which to evaluate the quantiles.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)}, otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#'
#' @return An integer vector of the same length as \code{p}, containing the quantiles corresponding to each probability in \code{p}.
#' @details The quantile function finds the smallest integer \eqn{k} such that the CDF at \eqn{k} is greater than or equal to the specified probability \eqn{p}.
#'
#' The function is vectorized, meaning that \code{p}, \code{lambda}, and \code{theta} can be vectors. If any of these parameters are single values (length 1), they will be used for all evaluations. All vector inputs must either be of the same length or have a length of 1.
#'
#' @examples
#' # Find the 25th, 50th, and 75th percentiles with lambda = 2 and theta = 1
#' qgdpois(c(0.25, 0.5, 0.75), lambda = 10, theta = 1)
#'
#' # Vectorized parameters: different values for theta
#' qgdpois(c(0.25, 0.5, 0.75), lambda = c(1, 5.3, 6.5), theta = 8)
#'
#' # Scalar values
#' qgdpois(0.99, lambda = 24, theta = 1)
#'
#' @export
qgdpois <- function(p, lambda, theta, lower.tail = TRUE, log.p = FALSE) {
  # Handle log.p transformation
  if (log.p) {
    p <- exp(p)
  }

  # Handle lower.tail transformation
  if (!lower.tail) {
    p <- 1 - p
  }

  # Validate probabilities
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("All probabilities 'p' must be between 0 and 1.")
  }

  # Get the lengths of each input
  lengths <- c(length(p), length(lambda), length(theta))

  # Check if all lengths are the same or if any length is 1
  if (!all(lengths == max(lengths) | lengths == 1)) {
    stop("All inputs must either be of the same length or length 1.")
  }

  # Determine the length to use for vectorization
  max_length <- max(lengths)

  # Repeat values if necessary to match the maximum length
  p <- if (length(p) == 1) rep(p, max_length) else p
  lambda <- if (length(lambda) == 1) rep(lambda, max_length) else lambda
  theta <- if (length(theta) == 1) rep(theta, max_length) else theta

  # Apply the scalar function element-wise
  sapply(seq_along(p), function(i) {
    current_p <- p[i]
    current_lambda <- lambda[i]
    current_theta <- theta[i]

    # Handle edge cases
    if (is.na(current_p)) {
      return(NA)
    }

    if (current_p == 0) {
      return(0)
    }

    if (current_p == 1) {
      # Assuming the distribution is defined for all k, return Inf to indicate no finite quantile
      return(Inf)
    }

    k <- 0
    while (pgdpois(k, current_lambda, current_theta) < current_p) {
      k <- k + 1
      # To prevent infinite loops in edge cases
      if (k > 1e6) {
        stop("Unable to find quantile. Please check the input parameters.")
      }
    }
    return(k)
  })
}











# R/logLikgd.R

#' Log-Likelihood Function for a sample of the GD-Poisson Distribution
#'
#' Calculates the log-likelihood for observed count data under the GD-Poisson distribution. This function is capable of handling both scalar and vector inputs for the observed counts (\code{k}), mean parameters (\code{lambda}), and dispersion parameters (\code{theta}). If any of these parameters are vectors, they must either be of the same length or have a length of 1, in which case the single value is used for all evaluations.
#'
#' @param k A vector of non-negative integers representing the observed count data.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#'
#' @return A single numeric value representing the sum of the log-likelihoods for the provided data and parameters.
#'
#' @details
#' The GD-Poisson (Poisson-Negative Zero) distribution is a generalization of the Poisson distribution that introduces arbitrary underdispersion. The log-likelihood is computed by first calculating the probability mass function (PMF) for each observed count \eqn{k} using the \code{dgdpois_scalar} function. To ensure numerical stability and avoid taking the logarithm of zero or \code{NA} values, any PMF values that are less than or equal to zero or \code{NA} are replaced with a very small positive number (\eqn{\epsilon = 1 \times 10^{-64}}).
#'
#'
#' @examples
#' # Single set of parameters
#' k_obs <- 3
#' lambda <- 2
#' theta <- 1
#' log_likelihood <- logLikgd(k_obs, lambda, theta)
#' print(log_likelihood)
#'
#' # Vector of observed counts with scalar lambda and theta
#' k_obs <- 0:5
#' log_likelihood <- logLikgd(k_obs, lambda = 2, theta = 1)
#' print(log_likelihood)
#'
#' # Vectors of the same length for k, lambda, and theta
#' k_obs <- 0:5
#' lambda <- c(2, 2, 3, 3, 4, 4)
#' theta <- c(1, 1.1, 1.2, 1.3, 1.4, 1.5)
#' log_likelihood <- logLikgd(k_obs, lambda, theta)
#' print(log_likelihood)
#'
#' # Example with invalid input lengths (should throw an error)
#' tryCatch({
#'   k_obs <- 0:5
#'   lambda <- c(2, 3)  # Length 2
#'   theta <- 1          # Length 1
#'   log_likelihood <- logLikgd(k_obs, lambda, theta)
#' }, error = function(e) {
#'   message("Error: ", e$message)
#' })
#'
#' @export
logLikgd <- function(k, lambda, theta) {
  # Get the lengths of each input
  lengths <- c(length(k), length(lambda), length(theta))

  # Check if all lengths are the same or if any length is 1
  if (!all(lengths == max(lengths) | lengths == 1)) {
    stop("All inputs must either be of the same length or length 1.")
  }
  # Determine the length to use for vectorization
  max_length <- max(lengths)

  # Repeat values if necessary to match the maximum length
  k <- if (length(k) == 1) rep(k, max_length) else k
  lambda <- if (length(lambda) == 1) rep(lambda, max_length) else lambda
  theta <- if (length(theta) == 1) rep(theta, max_length) else theta

  # Apply the function over the vectors
  Li=sapply(seq_along(k), function(i) {
    if (lambda[i] == 0) {
      1
    } else {
      dgdpois_scalar(k[i], lambda[i], theta[i])
    }
  })
  epsilon <- 1e-64
  Li[Li <= 0 | is.na(Li)] <- epsilon

  return(sum(log(Li)))
}


