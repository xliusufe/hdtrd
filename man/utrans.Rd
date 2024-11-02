\name{utrans}
\alias{utrans}
\title{
  Estimation of coefficient for the target data by transfer learning from the source data
}
\description{
  Provide the estimator of coefficient for the target data by transfer learning from the source data (Liu (2024)).
}

\usage{
  utrans(target, source, family = "gaussian", idtrans = NULL)
}

\arguments{
  \item{target}{The target dataset, a list, including \eqn{Y} (response), \eqn{\boldsymbol{X}} (covariates).}

  \item{source}{The source dataset, a list with sublist. Each sublist includes \eqn{Y} (response), \eqn{\boldsymbol{X}} (covariates). \code{source} could be \code{NULL}, in which case \code{utrans} only fits the target data by \code{glmnet}.}

  \item{family}{Family for generalized linear models, including `gaussian', `binomial', and `poisson'. Default is \code{family = "gaussian"}.}

  \item{idtrans}{The transferable source indices. It can be either a subset of {1,..., length(source)}. Default is \code{idtrans = NULL},  which is \code{idtrans = seq(length(source))}.}
}

\value{
  \item{fitglmnet}{The object from fitting \code{cv.glmnet} by CV method, see details in R package "glmnet".}

  \item{beta}{The coefficient (including intercept term) of the GLMs to fit target data by transfer learning.}

  \item{family}{The response type.}
}

\details{
See details in the paper Liu (2024)
}



\references{
Liu, S. (2024). Unified Transfer Learning Models for High-Dimensional Linear Regression. Proceedings of The 27th International Conference on Artificial Intelligence and Statistics, PMLR. 238, 1036-1044.

Liu, X. (2024). Testing relevant difference in high-dimensional linear regression with applications to detect transferability. Manuscript.

Tian, X., Lu, Y., and Li, W. (2015). A robust test for sphericity of high-dimensional covariance matrices. Journal of Multivariate Analysis, 141, 217-227.

Van den Meersche, K., Soetaert, K., and Van Oevelen, D. (2009). xsample(): An R Function for Sampling Linear Inverse Problems. Journal of Statistical Software, Code Snippets, 30, 1-15.
}

\examples{
  data(simulData_trans_gauss)
  fit <- utrans(target = dataset[[1]], source = dataset[-1], idtrans = seq(5))
  fit$beta[1:10]

}