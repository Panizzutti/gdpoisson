
logLikgd_scalar <- function(k, lambda, theta) {
  # Handle zero or negative PMF values to avoid log(0) or log of negative numbers
  if (lambda == 0) {
    0
  } else {
    p=dgdpois_scalar(k, lambda, theta)
    if (p <= 0 | is.na(p)) {return(-147)}
    log(p)
  }
}


neg_log_likelihood_loglink <- function(params, X, y) {
  # Extract beta coefficients and theta

  p <- length(params) - 1
  betas <- params[1:p]
  logtheta <- params[p+1]
  # Compute linear predictor and lambda
  eta <- X %*% betas
  lambda <- exp(eta)

  theta <-exp(logtheta)

  # Use sapply to compute dgdpois for each observation
  LLi <- sapply(seq_along(y), function(i) {
    logLikgd_scalar(y[i],lambda[i],theta)
  })

  negLL <- -sum(LLi)
  return(negLL)
}

MLE_find <- function(X, y, start_params = NULL, method = "L-BFGS-B", max_retries = 5) {
  n <- nrow(X)
  p <- ncol(X)

  # Generate initial starting values if not provided
  generate_starting_values <- function() {

    glm_fit <- glm(y ~ X-1, family = poisson(link = "log"))

    start_betas <- coef(glm_fit)
    pearson_residuals <- residuals(glm_fit, type = "pearson")
    start_logtheta <- log(sum(pearson_residuals^2) / (n - p))
    c(start_betas, start_logtheta)
  }

  # Get initial estimates from the Poisson model
  initial_estimates <- if (is.null(start_params)) {
    generate_starting_values()
  } else {
    start_params
  }

  if (length(initial_estimates) != p + 1) {
    stop("Length of start_params must be equal to number of betas plus one (theta).")
  }

  # Set lower and upper bounds
  #lower_bounds <- c(rep(-Inf, p), -Inf)  # Theta > 0, logtheta in R
  #upper_bounds <- c(rep(Inf, p), Inf)



  # Optimization with retry mechanism
  attempt <- 0
  while (attempt <= max_retries-1) {
    # Initialize var_cov_matrix and std_errors at the beginning of each attempt
    var_cov_matrix <- NULL
    std_errors <- rep(NA, p + 1)


    start_params <- initial_estimates + 0.1*attempt*rnorm(p+1)

    fit <- tryCatch({
      optim(
        par = start_params,
        fn = neg_log_likelihood_loglink,
        X = X,
        y = y,
        method = method,
        #lower = lower_bounds,
        #upper = upper_bounds,
        control = list(
          maxit = 1000,
          parscale = abs(start_params),
          pgtol = 1e-8,
          REPORT = 1  # Reports progress every iteration
        ),
        hessian = TRUE
      )
    }, error = function(e) {
      NULL
    })

    # Check if optimization was successful
    if (!is.null(fit) && fit$convergence == 0 && is.finite(fit$value)) {
      # Attempt to compute the variance-covariance matrix
      if (!is.null(fit$hessian)) {
        var_cov_matrix <- solve(fit$hessian)
        if (!is.null(var_cov_matrix)) {
          # Hessian inversion successful
          std_errors <- sqrt(diag(var_cov_matrix))
          # Optimization and Hessian inversion successful, exit loop
          break
        } else {
          # Singular Hessian, retry optimization
          if (attempt == max_retries - 1) warning("Hessian is singular. Final attempt failed.")
        }
      } else {
        # Hessian not available, retry optimization
        if (attempt == max_retries - 1) warning("Hessian not available. Final attempt failed.")
      }
    } else {
      # Optimization failed, retry with new starting values
      if (attempt == max_retries - 1) warning("Optimization failed after final attempt. Review data or model specification.")
    }

    attempt <- attempt + 1
  }


  # If optimization still failed after max_retries
  if (attempt > max_retries || is.null(fit) || fit$convergence != 0 || !is.finite(fit$value) || is.null(var_cov_matrix)) {
    stop("Optimization failed after multiple attempts. Consider checking the data or model specification.")
  }

  fit$std_errors= std_errors
  fit$var_cov_matrix= var_cov_matrix

  return(fit)

}

mle_gdpois_loglink <- function(X, y, start_params = NULL, method = "L-BFGS-B", max_retries = 5) {
  n <- nrow(X)
  p <- ncol(X)

  fit= MLE_find(X, y,start_params=start_params, method = method, max_retries = max_retries)

  NULLfit= MLE_find( as.matrix(rep(1,n),ncol=1), y,start_params=start_params, method = method, max_retries = max_retries)

  # Extract estimated parameters
  est_params <- fit$par
  betas_est <- est_params[1:p]
  theta_est <- exp(est_params[p + 1])

  logLikelihood = -fit$value

  #satLLi <- sapply(seq_along(y), function(i) {logLikgd_scalar(y[i],y[i],theta_est)})
  saturatedLL=0 #sum(satLLi)

  resdeviance= 2*(saturatedLL - logLikelihood)
  nullLL=-NULLfit$value
  #meannull=mean(y)
  #nullLLi <- sapply(seq_along(y), function(i) {
    #logLikgd_scalar(y[i],meannull,theta_est)
  #})
  #nullLL= sum(nullLLi)

  nulldeviance=2*(saturatedLL-nullLL)

  # Compile results
  result <- list(
    coefficients = betas_est,
    theta = theta_est,
    std_errors = fit$std_errors,
    logLik = logLikelihood,
    aic = 2 * (p+1) - 2 * logLikelihood,
    bic = log(n) * (p+1) - 2 * logLikelihood,
    null.deviance = nulldeviance,
    deviance = resdeviance,
    df.null=n-2,
    df.residual=n-p-1,
    converged = fit$convergence,
    var_cov_matrix = fit$var_cov_matrix,
    method= method,
    X = X,
    y = y,
    n = n,
    p = p+1,
    fitted.values = exp(X%*%betas_est),
    residuals = y-exp(X%*%betas_est),
    iter = fit$counts["function"]
  )

  return(result)
}


# R/glm_gdp.R

#' Fit a GD-Poisson Generalized Linear Model using a logarithmic link function.
#'
#' Fits a generalized linear model using the GD-Poisson distribution for the response variable.
#'
#' @param formula An object of class \code{\link{formula}}. The model formula specifying the response and predictors.
#' @param data A data frame containing the variables specified in the formula.
#' @param method Character string specifying the optimization method to be used in the function. Default is \code{"L-BFGS-B"}.
#' @param max_retries Integer specifying the maximum number of retries for the optimization algorithm.
#' @return An object of class \code{glm.gdp} containing the fitted model results, including coefficients, fitted values, and other relevant information.
#'
#' @details
#' The \code{glm.gdp} function fits a GD-Poisson generalized linear model to the provided data. By specifying \code{method} users can control the optimization process used in maximum likelihood estimation.
#'
#' @examples
#' # Fit a GD-Poisson model with default settings
#' fit <- glm.gdp(y ~ x1 + x2, data = my_data)
#'
#' # Fit a GD-Poisson model with custom optimization method and maximum retries
#' fit <- glm.gdp(y ~ x1 + x2, data = my_data, method = "BFGS", max_retries = 10)
#'
#' @export
glm.gdp <- function(formula, data, method = "L-BFGS-B", max_retries = 5) {
  # Extract response and design matrix
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)
  start_params = NULL
  # Fit the model using mle_gdppois_loglink with specified parameters
  fit_result <- mle_gdpois_loglink(X, y, start_params = start_params, method = method, max_retries = max_retries)

  # Store additional information
  fit_result$call <- match.call()
  fit_result$formula <- formula
  fit_result$data <- data

  # Set the class of the object
  class(fit_result) <- "glm.gdp"

  return(fit_result)
}



#' Summary Method for GD-Poisson GLM Objects
#'
#' Provides a comprehensive summary of a fitted GD-Poisson generalized linear model, including coefficients,
#' standard errors, z-values, p-values, residual summaries, dispersion parameter, log-likelihood,
#' AIC, BIC, deviance measures, degrees of freedom, number of iterations, and convergence status.
#'
#' @param object An object of class \code{gdpois_glm} produced by \code{\link{glm.gdp}}.
#' @param ... Additional arguments (currently not used).
#'
#' @return An object of class \code{summary.gdp_glm} containing the summary information.
#'
#' @export
summary.glm.gdp <- function(object, ...) {

  # Create a coefficients table
  coefficients_table <- data.frame(
    Estimate     = object$coefficients,
    `Std. Error` = object$std_errors[1:length(object$coefficients)],
    `z value`    = object$coefficients / object$std_errors[1:length(object$coefficients)],
    `Pr(>|z|)`   = 2 * (1 - pnorm(abs(object$coefficients / object$std_errors[1:length(object$coefficients)]))),
    check.names  = FALSE
  )

  # Assign parameter names as row names
  rownames(coefficients_table) <- colnames(object$X)

  # Compile the summary as a list directly referencing object components
  summary_output <- list(
    call            = object$call,
    terms           = object$formula,
    coefficients     = coefficients_table,
    dispersion       = object$theta,
    null.deviance    = object$null.deviance,
    residual.deviance= object$deviance,
    aic              = object$aic,
    bic              = object$bic,
    df.null          = object$df.null,
    df.residual      = object$df.residual,
    logLik           = object$logLik,
    iter             = object$iter,
    fitted.values    = object$fitted.values,
    data             = object$data,
    n                = object$n,
    p                = object$p,  # Including theta
    converged        = object$converged
  )

  # Assign class to the summary object
  class(summary_output) <- "summary.glm.gdp"

  return(summary_output)
}




#' Print Method for Summary GD-Poisson GLM Objects
#'
#' Prints a detailed summary of a fitted GD-Poisson generalized linear model to the console,
#' including call, residual summaries, coefficients table, dispersion parameter, deviance measures,
#' information criteria, degrees of freedom, and number of iterations.
#'
#' @param x An object of class \code{summary.glm.gdp}.
#' @param ... Additional arguments (currently not used).
#'
#' @return Prints the summary to the console.
#'
#' @export
print.summary.glm.gdp <- function(x, ...) {


  # Print the Call
  cat("\nCall:\n")
  print(x$call)

  # Print Coefficients
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = 4, signif.stars = TRUE)

  # Print Dispersion Parameter
  cat("\nDispersion parameter (theta):", formatC(x$dispersion, digits = 4), "\n")

  # Print Null and Residual Deviance
  cat("\nNull deviance:", formatC(x$null.deviance, digits = 4),
      "on", x$df.null, "degrees of freedom\n")
  cat("Residual deviance:", formatC(x$residual.deviance, digits = 4),
      "on", x$df.residual, "degrees of freedom\n")

  # Print AIC and BIC
  cat("AIC:", formatC(x$aic, digits = 4), "\n")

  # Print Number of Iterations
  cat("Number of Fisher Scoring iterations:", x$iter, "\n")

  # Check for Convergence Status
  if (!is.null(x$converged)) {
    if (x$converged == 0) {
      cat("Optimization converged successfully.\n")
    } else {
      cat("Optimization did not converge.\n")
    }
  } else {
    cat("Convergence information not available.\n")
  }
}



