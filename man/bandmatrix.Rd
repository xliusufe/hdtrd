\name{bandmatrix}
\alias{bandmatrix}
\title{
  Construct a sparse banded matrix.
}
\description{
  Provide a sparse banded matrix.
}

\usage{
  bandmatrix(rho, p, T = 5)
}

\arguments{
  \item{rho}{A vector with length \eqn{T}.}

  \item{p}{The dimension of the banded matrix.}

  \item{T}{The width of band. Default is \code{T = 5}.}
}

\value{
  \item{sighalf}{The matrix \eqn{\Gamma \in \mathcal{R}^{(p+T)\times p}} satisfying \eqn{\Sigma = \Gamma^T\Gamma}.}

  \item{sigma}{The sparse banded matrix \eqn{\Sigma \in \mathcal{R}^{p\times p}}.}
}

\references{
Chen, Z., Cheng, V. X. and Liu, X. (2024). Hypothesis testing on high dimensional quantile regression. Journal of Econometrics.
}

\examples{
  p   <- 6
  T   <- 3
  rho <- seq(T)/(T+1)
  fit <- bandmatrix(rho, p, T)
  fit$sigma
}