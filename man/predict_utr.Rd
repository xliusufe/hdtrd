\name{predict_utr}
\alias{predict_utr}
\title{
  Prediction of a new predictor
}
\description{
  Provide the prediction for a new predictor.
}

\usage{
  predict_utr(fittrans, X, type = "response")
}

\arguments{
  \item{fittrans}{An object from fitting \code{utrans}.}

  \item{X}{A new predictor, a matrix with dimension \eqn{n\times p}.}

  \item{type}{The type of prediction, including "response" (Default) and "class".  Here "response" provides the predicted probability when \code{family = "binomial"}. "class" predict \code{0/1} response for logistic regression. Applies only when \code{family = "binomial"}.}
}

\value{
  \item{yhat}{The new response \eqn{\hat{y}} based on the new predictor \eqn{x}.}
}

\details{
See details in the paper Liu (2024).
}



\references{
Liu, S. (2024). Unified Transfer Learning Models for High-Dimensional Linear Regression. Proceedings of The 27th International Conference on Artificial Intelligence and Statistics, PMLR. 238, 1036-1044.

Liu, X. (2024). Testing relevant difference in high-dimensional linear regression with applications to detect transferability. Manuscript.

Tian, X., Lu, Y., and Li, W. (2015). A robust test for sphericity of high-dimensional covariance matrices. Journal of Multivariate Analysis, 141, 217-227.

Van den Meersche, K., Soetaert, K., and Van Oevelen, D. (2009). xsample(): An R Function for Sampling Linear Inverse Problems. Journal of Statistical Software, Code Snippets, 30, 1-15.
}

\examples{
  data(simulData_trans_gauss)
  fittrans <- utrans(target = dataset[[1]], source = dataset[-1], idtrans = seq(5))

  p = ncol(dataset[[1]]$X)
  n = 5
  sig = toeplitz(0.5^(c(1:p)-1))
  sighalf = chol( sig )
  x = matrix(rnorm(n*p), nrow = n) %*% sighalf

  predict_utr(fittrans, x)
}