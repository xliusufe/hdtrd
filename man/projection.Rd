\name{projection}
\alias{projection}
\title{
  Projection of \eqn{y} onto the closure of covariates \eqn{x}
}
\description{
  Provide the projection of \eqn{y} onto the closure of covariates \eqn{x}.
}

\usage{
  projection(x, y, family = "gaussian", method = 'lasso', isresid = TRUE)
}

\arguments{
  \item{x}{Covariates, a \eqn{n\times p}-matrix.}

  \item{y}{Response, a \eqn{n}-vector.}

  \item{family}{Family for the generalized linear models, including `gaussian', `binomial', and `poisson'. Default is \code{family = "gaussian"}.}

  \item{method}{There are two methods, "qfabs" and "lasso", to estimate the nuisance parameter \eqn{\boldsymbol{\alpha}} in quantile regression. Default is \code{method = 'lasso'}.}

  \item{isresid}{logical. Projected residual \eqn{\hat{\boldsymbol{\eta}}=\boldsymbol{x} - \hat{H}\boldsymbol{z}} is output if \code{isresid = TRUE}. Coefficient matrix \eqn{\hat{H}} is calculated if \code{isresid = FALSE}. Default is \code{resids = TRUE}.}
}

\value{
  \item{proj}{Projection.}
}

\details{
  High-dimensional test of relevant difference and its application to transferability test in the generalized Linear regression models (see details in the paper Liu (2024))
  \deqn{
  y_i = H\boldsymbol{X}_{i}^T.
  }
}


\references{
Cheng, C., Feng, X., Huang, J. and Liu, X. (2022). Regularized projection score estimation of treatment effects in high-dimensional quantile regression. Statistica Sinica. 32, 23-41.

Liu, X. (2024). Testing relevant difference in high-dimensional linear regression with applications to detect transferability. Manuscript.
}

\examples{
  data(simulData_test_gauss)
  x <- datahb$X
  y <- datahb$Y
  proj <- projection(x, y)
}