% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdp_regression.R
\name{summary.glm.gdp}
\alias{summary.glm.gdp}
\title{Summary Method for GD-Poisson GLM Objects}
\usage{
\method{summary}{glm.gdp}(object, ...)
}
\arguments{
\item{object}{An object of class \code{gdpois_glm} produced by \code{\link{glm.gdp}}.}

\item{...}{Additional arguments (currently not used).}
}
\value{
An object of class \code{summary.gdp_glm} containing the summary information.
}
\description{
Provides a comprehensive summary of a fitted GD-Poisson generalized linear model, including coefficients,
standard errors, z-values, p-values, residual summaries, dispersion parameter, log-likelihood,
AIC, BIC, deviance measures, degrees of freedom, number of iterations, and convergence status.
}
