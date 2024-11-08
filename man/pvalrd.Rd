\name{pvalrd}
\alias{pvalrd}
\title{
  P-value for high-dimensional testing of relevant difference in the generalized linear regression models when the nuisance parameter is high-dimensional
}
\description{
  Provide p-value for high-dimensional testing of relevant difference in generalized linear regression models (Liu (2024)) when the nuisance parameter is high-dimensional.
}

\usage{
  pvalrd(data, family = "gaussian", delta0 = 0.1, method = 'lasso',
            resids = NULL, sigma2 = NULL, lammax = NULL)
}

\arguments{
  \item{data}{A list, including \eqn{Y} (response), \eqn{\boldsymbol{X}}, \eqn{\boldsymbol{Z}}, where \eqn{\boldsymbol{Z}} is high-dimensional.}

  \item{family}{Family for the generalized linear models, including `gaussian', `binomial', and `poisson'. Default is \code{family = "gaussian"}.}

  \item{delta0}{Relevant difference, a given value by hypothesis test problem \eqn{H_0: \|\beta\|_2\leq\delta_0}. Default is \code{delta0 = 0.1}.}

  \item{method}{There are two methods, "qfabs" and "lasso", to estimate the nuisance parameter \eqn{\boldsymbol{\alpha}} in quantile regression. Default is \code{method = 'lasso'}.}

  \item{resids}{An \eqn{n}-vector, which is residual of the GLM under \eqn{H_0}. Default is \code{resids = NULL}, where the canonical link function is used if \code{resids} and \code{psi} are \code{NULL}.}

  \item{sigma2}{Estimator of error's variance if family = "gaussian". Default is \code{sigma2 = NULL}, where \code{sigma2 = 1}. }

  \item{lammax}{Estimator of the largest eigenvalue \eqn{\sup_{\|\beta\|_2\leq\delta_0}\beta^T\Sigma^2\beta}. Default is \code{lammax = NULL}, which is estimated empirically by \eqn{\lambda_{\max}(S_n)/(1+\sqrt{p/n})}, see details in Liu (2024). }
}

\value{
  \item{pvals}{P-value of the corresponding test statistic.}

  \item{Tn}{Standardized test statistic.}
}

\details{
  High-dimensional test of relevant difference and its application to transferability test in the generalized Linear regression models (see details in the paper Liu (2024))
  \deqn{
  \mu_i = \boldsymbol{X}_{i}^T\boldsymbol{\beta} +\boldsymbol{Z}_i^T\boldsymbol{\gamma},
  }
where \eqn{\boldsymbol{Z}^T\boldsymbol{\gamma}} is the control mean function, and \eqn{\boldsymbol{X}} is high-dimensional.

The hypothesis test problem is
\deqn{
  H_0: \|\boldsymbol{\beta}\|\leq \delta_0 \quad versus\quad H_1: \|\boldsymbol{\beta}\|> \delta_0.
}
}


\references{
Karoui, N, E. (2008) Spectrum estimation for large dimensional covariance matrices using random matrix theory. The Annals of Statistics, 36(6), 2757-2790.

Kong, W. and Valiant, G. (2017). Spectrum estimation from samples. Annals of Statistics. 45, 2218-2247.

Liu, X. (2024). Testing relevant difference in high-dimensional linear regression with applications to detect transferability. Manuscript.

Tian, X., Lu, Y., and Li, W. (2015). A robust test for sphericity of high-dimensional covariance matrices. Journal of Multivariate Analysis, 141, 217-227.

Van den Meersche, K., Soetaert, K., and Van Oevelen, D. (2009). xsample(): An R Function for Sampling Linear Inverse Problems. Journal of Statistical Software, Code Snippets, 30, 1-15.
}

\examples{
  data(simulData_test_gauss)
  pvals <- pvalrd(data = datahb)
  pvals$pvals

}