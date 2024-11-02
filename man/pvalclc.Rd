\name{pvalclc}
\alias{pvalclc}
\title{
  P-value for high-dimensional test in the generalized linear regression models
}
\description{
  Provide p-value for high-dimensional test in the generalized linear regression models when the nuisance parameter is high-dimensional, see Chen et. al. (2022) for details.
}

\usage{
  pvalclc(data, family = 'gaussian', method = 'lasso', resids = NULL, psi = NULL)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{\boldsymbol{X}}, \eqn{\boldsymbol{Z}}, where \eqn{\boldsymbol{X}} can be \code{NULL}.}

  \item{family}{Family for the generalized linear models, including `gaussian', `binomial', and `poisson'. Default is \code{family = "gaussian"}.}

  \item{method}{There are two methods, "gfabs" and "lasso", to estimate the nuisance parameter \eqn{\boldsymbol{\alpha}} in GLMs. Default is \code{method = 'lasso'}, which calls \code{glmnet}.}

  \item{resids}{An \eqn{n}-vector, which is residual of the GLM. Default is \code{resids = NULL}. The canonical link function is used if \code{resids} and \code{psi} are \code{NULL}.}

  \item{psi}{An \eqn{n}-vector, which is \eqn{\psi(X_i,\beta_0,\phi)=g'(X_i^{\top}\beta_0)/V(\mu_i(\beta_0);\phi)}, see Guo and Chen (2016) for the details. Default is \code{psi = NULL}.
    The canonical link function is used if both \code{resids} and \code{psi} are \code{NULL}. Here, \code{psi = rep(1,n)} if \code{psi = NULL}.}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}

  \item{Tn}{Test statistic \eqn{|\hat{U}_n|/\sqrt{2\hat{R}_n}}. Reject \eqn{H_0} if \eqn{|\hat{U}_n|/\sqrt{2\hat{R}_n}>z_{1-\alpha/2}}.}
}

\details{
  The generalized linear regression models (see details in the paper Guo and Chen (2016))
  \deqn{
  \mu_i = \boldsymbol{X}_{i}^T\boldsymbol{\alpha} +\boldsymbol{Z}_i^T\boldsymbol{\beta},
  }
where \eqn{\boldsymbol{X}^T\boldsymbol{\alpha}} is a baseline mean function.

The hypothesis test problem is
\deqn{
  H_0: \boldsymbol{\beta}=\boldsymbol{0}\quad versus\quad H_1: \boldsymbol{\beta}\neq \boldsymbol{0}.
}

One can input estimated residual \code{resids = } \eqn{y_i - \hat{\mu}_i} and \code{psi = NULL} which produces the test statistic and p-value given by Chen et. al. (2022),
  where \eqn{\hat{\mu}_i} is an estimator of \eqn{\mu_i} according to Chen et. al. (2023).
}



\references{

Guo, B.and Chen, S. X. (2016). Tests for high dimensional generalized linear models. Journal of the Royal Statistical Society, Series B, 78, 1079-1102.

Chen, J., Li, Q., and Chen, H. Y. (2023). Testing generalized linear models with highdimensional nuisance parameters. Biometrika. 110, 83-99.
}

\examples{
  data(simulData_test_gauss)
  pvals <- pvalclc(data = datahb, family = "gaussian")
  pvals$pvals
}