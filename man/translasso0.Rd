\name{translasso0}
\alias{translasso0}
\title{
  Estimation of coefficient for the target data by oracle Trans-lasso from the source data
}
\description{
  Provide the estimator of coefficient for the target data by oracle Trans-Lasso from the source data (Li et al. (2022b)).
}

\usage{
  translasso0(target, source = NULL, idtrans = NULL, nvec = NULL, lamconst = NULL, l1 = TRUE)
}

\arguments{
  \item{target}{The target dataset, a list, including \eqn{Y} (response), \eqn{\boldsymbol{X}} (covariates).}

  \item{source}{The source dataset, a list with sublist. Each sublist includes \eqn{Y} (response), \eqn{\boldsymbol{X}} (covariates). \code{source} could be \code{NULL}, in which case \code{utrans} only fits the target data by \code{glmnet}.}

  \item{idtrans}{The transferable source indices. It can be either a subset of {1,..., length(source)}. Default is \code{idtrans = NULL},  which is \code{idtrans = seq(length(source))}.}

  \item{nvec}{A vector integers with length \eqn{K_s+1}, each element of which is the number of indices of samples. Default is \code{idtrans = NULL},  which is \code{nvec = c(n0,n1,...,n_K)}.}

  \item{lamconst}{A numeric number, which is a constant in the form \code{lamconst*sqrt(2*log(p)/n)}. Default is \code{lamconst = NULL}, which means that \code{lamconst} is selected by \code{glmnet}.}

  \item{l1}{method to estimate the parameter in linear regression. Default is \code{l1 = TRUE}, which means that lasso is applied.}
}

\value{
  \item{beta}{The coefficient (including intercept term) fitted target data by trans-lasso.}
}

\details{
See details in the paper Li et al. (2022b)
}



\references{
Li, S., Cai, T. T., and Li, H. (2022a). Transfer Learning in Large-Scale Gaussian Graphical Models with False Discovery Rate Control. Journal of the American Statistical Association, 118, 2171-2183.

Li, S., Cai, T. T., and Li, H. (2022b). Transfer Learning for High-Dimensional Linear Regression: Prediction, Estimation and Minimax Optimality. Journal of the Royal Statistical Society Series B, 84, 149-173.

Liu, S. (2024). Unified Transfer Learning Models for High-Dimensional Linear Regression. Proceedings of The 27th International Conference on Artificial Intelligence and Statistics, PMLR. 238, 1036-1044.

Liu, X. (2024). Testing relevant difference in high-dimensional linear regression with applications to detect transferability. Manuscript.

Tian, X., Lu, Y., and Li, W. (2015). A robust test for sphericity of high-dimensional covariance matrices. Journal of Multivariate Analysis, 141, 217-227.

Van den Meersche, K., Soetaert, K., and Van Oevelen, D. (2009). xsample(): An R Function for Sampling Linear Inverse Problems. Journal of Statistical Software, Code Snippets, 30, 1-15.
}

\examples{
  data(simulData_trans_gauss)
  fit <- translasso0(target = dataset[[1]], source = dataset[-1], idtrans = seq(5))
  fit$beta[1:10]

}