\name{eigmax}
\alias{eigmax}
\title{
  Estimation of the largest eigenvalue of covariance of a high-dimensional vector
}
\description{
  Provide the estimator of the largest eigenvalue of covariance of a high-dimensional vector (Liu (2024)), as well as all estimated eigenvalues.
}

\usage{
  eigmax(X, zK = NULL, tJ = NULL, K = 1000, J = 1000, method = 'mpmo',
            nmoms = NULL, timeout = 0L)
}

\arguments{
  \item{X}{A data matrix with dimension \eqn{n\times p}.}

  \item{zK}{A matrix with dimension \eqn{K\times 2}, a given complex number, where the first column is the real part and the second column is the imaginary part. Default is \code{zK = NULL}, where \code{zK[,1] = rnorm(K)} is generated from standard normal distribution, and \code{zK[,2] =  rep(1,K)/sqrt(n)}.}

  \item{tJ}{A \eqn{J}-vector. Default is \code{tJ = NULL}, where \code{tJ} is a grid of points in the interval \eqn{[\lambda_{\min}(\Sigma),\lambda_{\max}(\Sigma)]}.}

  \item{K}{A positive integer, which is the number of complex numbers \code{zK}. Default is \code{K = 1000}.}

  \item{J}{A positive integer, which is the length of \code{tJ}. Default is \code{J = 1000}.}

  \item{method}{There are three methods, \code{'mpmo'}, \code{'mplp'} and \code{'empi'}, to estimate the largest eigenvalue of \eqn{\Sigma}, see details in Liu (2024). Default is \code{method = 'mpmo'}.}

  \item{nmoms}{The number of moments.  Default is \code{nmoms = NULL}, where \code{nmoms = 7} if \code{method = 'mpmo'}, \code{nmoms = 4} if \code{method = 'mplp'}, and \code{nmoms} is useless if \code{method = 'empi'}. }

  \item{timeout}{An integer: timeout variable in seconds, defaults to 0L which means no limit is set, see details in the function \code{linp} of R package "limSolve".}
}

\value{
  \item{lammax}{Estimator of the largest eigenvalue of \eqn{\Sigma}.}

  \item{lamest}{All estimated eigenvalues of \eqn{\Sigma}.}
}

\details{
See details in the paper Liu (2024).

Here, for the methods to estimate the largest eigenvalue of \eqn{\Sigma}, \code{'mpmo'} denotes the MPMO method; \code{'mplp'} denotes the MPLP method; and \code{'empi'} denotes the EMPI method.
}



\references{
Karoui, N, E. (2008) Spectrum estimation for large dimensional covariance matrices using random matrix theory. The Annals of Statistics, 36(6), 2757-2790.

Kong, W. and Valiant, G. (2017). Spectrum estimation from samples. Annals of Statistics. 45, 2218-2247.

Liu, X. (2024). High-dimensional test of relevant difference and its application totransfer learning. Manuscript.

Tian, X., Lu, Y., and Li, W. (2015). A robust test for sphericity of high-dimensional covariance matrices. Journal of Multivariate Analysis, 141, 217-227.

Van den Meersche, K., Soetaert, K., and Van Oevelen, D. (2009). xsample(): An R Function for Sampling Linear Inverse Problems. Journal of Statistical Software, Code Snippets, 30, 1-15.
}

\examples{
  p = 300
  n = 200
  sig = toeplitz(0.5^(c(1:p)-1))
  sighalf = chol( sig )
  X 	= matrix(rnorm(n*p), nrow = n) %*% sighalf
  eigens = eigmax(X = X, method = 'mpmo')
  eigens$lammax

}