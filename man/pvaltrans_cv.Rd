\name{pvaltrans_cv}
\alias{pvaltrans_cv}
\title{
  P-value for high-dimensional testing of relevant difference in high-dimensional transfer learning in the generalized linear regression models via cross validation method.
}
\description{
  Provide p-value for high-dimensional testing of relevant difference in high-dimensional transfer learning in the generalized linear regression models via cross validation method (Liu (2024)).
}

\usage{
  pvaltrans_cv(target, source, family = "gaussian", delta0 = 0.1, nsource = 10,
                method = 'lasso', ncv = 10, alpha = 0.05, resids = NULL,
                isproj = FALSE, proj = NULL,sigma2 = NULL, lammax = NULL,
                nmoms = NULL, zK = NULL, J = NULL, K = NULL, timeout = 0)
}


\arguments{
  \item{target}{The target dataset, a list, including \eqn{Y} (response), \eqn{\boldsymbol{X}} (covariates).}

  \item{source}{The source dataset, a list with sublist. Each sublist includes \eqn{Y} (response), \eqn{\boldsymbol{X}} (covariates).}

  \item{family}{Family for generalized linear models, including `gaussian', `binomial', and `poisson'. Default is \code{family = "gaussian"}.}

  \item{delta0}{Relevant difference, a given value by hypothesis test problem \eqn{H_0: \|\beta\|_2\leq\delta_0}. Default is \code{delta0 = 0.1}.}

  \item{nsource}{The number of source datasets. Default is \code{nsource = 10}.}

  \item{method}{There are two methods, "glm" and "lasso", to estimate the nuisance parameter \eqn{\boldsymbol{\alpha}} under the null hypothesis in  the generalized linear regression models, where "glm" method estimates nuisance parameter for classic low-dimensional setting, and "lasso" for high-dimensional setting. Default is \code{method = 'lasso'} for high-dimensional setting.}

  \item{ncv}{Number of folds in the cross-validation, which is used to select transferable level \eqn{\delta_0}. Default is \code{ncv = 10}.}

  \item{alpha}{Significant level, which is used to select transferable level \eqn{\delta_0}. Default is \code{alpha = 0.05}.}

  \item{resids}{An \eqn{n}-vector, which is residual of GLM under \eqn{H_0}. Default is \code{resids = NULL}, where the canonical link function is used if \code{resids} and \code{psi} are \code{NULL}.}

  \item{isproj}{logical. Projection score method is applied if \code{isproj = TRUE}. Default is \code{isproj = FALSE}, which means that no projection score is applied.}

  \item{proj}{The estimated residual of projection score, a list, where each element is a \eqn{n\times p}-matrix, \eqn{\hat{\eta}=x-\hat{H}z}. Default is \code{proj = NULL}, which means that projection score is calculated. }

  \item{sigma2}{Estimator of error's variance if family = "gaussian". Default is \code{sigma2 = NULL}, where \code{sigma2 = 1}. }

  \item{lammax}{Esimator of the largest eigenvalue \eqn{\sup_{\|\beta\|_2\leq\delta_0}\beta^T\Sigma^2\beta}, see details in \code{eigmax}. Default is \code{lammax = NULL}, where \code{lammax} is estimated by EMPI method, see \code{eigmax}. If \code{testmethd = 'pvalrd'}, there are two choices \code{lammax = 'mpmo'} or  \code{lammax = 'mplp'}. It is useless if \code{testmethd = 'pvalclc'}. }

  \item{nmoms}{The number of moments.  Default is \code{nmoms = NULL}, where \code{nmoms = 7} if \code{method = 'mpmo'}, \code{nmoms = 4} if \code{method = 'mplp'}, and \code{nmoms} is useless if \code{method = 'empi'}. }

  \item{zK}{A matrix with dimension \eqn{K\times 2}, a given complex number, where the first column is the real part and the second column is the imaginary part. Default is \code{zK = NULL}, where \code{zK[,1] = rnorm(K)} is generated from standard normal distribution, and \code{zK[,2] =  rep(1,K)/sqrt(n)}.}

  \item{J}{A positive integer, which is the length of \code{tJ}. Default is \code{J = NULL}, which means \code{J = max(500,3*n,2*p)+200}.}

  \item{K}{A positive integer, which is the number of complex numbers \code{zK}. Default is \code{K = NULL}, which means \code{K = max(500,3*n,2*p)+200}.}

  \item{timeout}{An integer: timeout variable in seconds, defaults to 0L which means no limit is set, see details in the function \code{linp} of R package "limSolve".}
}

\value{
  \item{pvals}{P-value of the corresponding test statistic, which is a vector with length \code{nsource}.}
  \item{s_opt}{The \code{s_opt}th \eqn{\delta_0} is Selected.}
}

\details{
  High-dimensional test of relevant difference and its application to transferability test in the generalized Linear regression models (see details in the paper Liu (2024)).

  Linear regression model for target data:
  \deqn{
  Y_{0i} = \boldsymbol{X}_{0i}^T\boldsymbol{\beta}_0 + \epsilon_{0i},
  }

and 

  linear regression model for the \eqn{k}th source data:
  \deqn{
  Y_{ki} = \boldsymbol{X}_{ki}^T\boldsymbol{\beta}_k + \epsilon_{ki},
  }
where \eqn{\boldsymbol{X}^T\boldsymbol{\beta}} is a baseline mean function, and \eqn{\boldsymbol{X}} is high-dimensional.

The hypothesis test problem is
\deqn{
  H_0: \|\boldsymbol{\beta}-\boldsymbol{\beta}_0\|\leq \delta_0 \quad versus\quad H_1: \|\boldsymbol{\beta}-\boldsymbol{\beta}_0\|> \delta_0.
}

Here, for the methods to estimate the largest eigenvalue of \eqn{\Sigma}, \code{'mpmo'} denotes the MPMO method; \code{'mplp'} denotes the MPLP method; and \code{'empi'} denotes the EMPI method.
}



\references{
Chen, Z., Cheng, V. X. and Liu, X. (2024). Hypothesis testing on high dimensional quantile regression. Journal of Econometrics.

Karoui, N, E. (2008) Spectrum estimation for large dimensional covariance matrices using random matrix theory. The Annals of Statistics, 36(6), 2757-2790.

Kong, W. and Valiant, G. (2017). Spectrum estimation from samples. Annals of Statistics. 45, 2218-2247.

Liu, S. (2024). Unified Transfer Learning Models for High-Dimensional Linear Regression. Proceedings of The 27th International Conference on Artificial Intelligence and Statistics, PMLR. 238, 1036-1044.

Liu, X. (2024). Testing relevant difference in high-dimensional linear regression with applications to detect transferability. Manuscript.

Liu, X., Zheng, S. and Feng, X. (2020). Estimation of error variance via ridge regression. Biometrika. 107, 481-488.

Tian, Y. and Feng, Y. (2023) Transfer Learning Under High-Dimensional Generalized Linear Models. Journal of the American Statistical Association, 118, 2684-2697.

Yang, W., Guo, X. and Zhu, L. (2023). Score function-based tests for ultrahigh-dimensional linear models. arXiv:2212.08446.

Zhang, X. and Cheng, G. (2017). Simultaneous inference for high-dimensional linear models. Journal of the American Statistical Association, 112, 757-768.

Van den Meersche, K., Soetaert, K., and Van Oevelen, D. (2009). xsample(): An R Function for Sampling Linear Inverse Problems. Journal of Statistical Software, Code Snippets, 30, 1-15.
}

\examples{
  data(simulData_trans_gauss)
  np      <- dim(dataset[[1]]$X)
  delta0 	<- c(1:10)*log(np[1])/np[2]
  ## pvals <- pvaltrans_cv(target = dataset[[1]], source = dataset[-1], delta0 = delta0, nsource = 1)
}