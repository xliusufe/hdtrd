\name{simulData}
\alias{simulData}
\docType{data}
\title{
Simulated data for generalized linear regression models
}
\description{
Simulated data for generalized linear regression models.
\itemize{
\item `Linear regression' for testing relevant difference (\code{simulData_test_gauss}),
\item `Poisson regression' for testing relevant difference (\code{simulData_test_poiss}), 
\item `Logistic regression' for testing relevant difference (\code{simulData_test_binom}).
\item `Linear regression' for transfer learning (\code{simulData_trans_gauss}),
\item `Poisson regression' for transfer learning (\code{simulData_trans_poiss}), and
\item `Logistic regression' for transfer learning (\code{simulData_trans_binom}).
}

Each dataset includes a list entitled 
\itemize{
\item \code{datahb} in \code{simulData_test} for linear regression models,
\item \code{data_binom} in \code{simulatedData_Binom} for logistic regression models,
\item \code{data_poiss} in \code{simulatedData_Poiss} for Poisson regression models, 
\item \code{dataset} in \code{simulData_trans_gauss}, \code{simulData_trans_binom} and \code{simulData_trans_poiss} for linear regression, logistic regression and Poisson regression models, respectively. \code{dataset[[1]]} is the target dataset, and \code{dataset[-1]} is the 10 source datasets.
}
}
\usage{
data(simulData_test_gauss)
}
\details{
For \code{simulData_test_gauss}, we simulated data generated from linear regression models
  \deqn{
  Y_i = \boldsymbol{X}_{i}^T\boldsymbol{\alpha} +\boldsymbol{Z}_i^T\boldsymbol{\beta} +\epsilon_i,
  }
where \eqn{\boldsymbol{X}^T\boldsymbol{\alpha}} is a baseline mean function.

    \itemize{
        \item Y: the response, an \eqn{n}-vector,

        \item X: the baseline variable with dimension \eqn{n\times p},

        \item Z: the interested variable with dimension \eqn{n\times q}.
    }

For \code{simulData_trans_gauss}, we simulated data generated from linear regression models
  \deqn{
  Y_{0i} = \boldsymbol{X}_{0i}^T\boldsymbol{\beta}_0 + \epsilon_{0i},
  }

and 

  Linear regression model for the \eqn{k}th source data:
  \deqn{
  Y_{ki} = \boldsymbol{X}_{ki}^T\boldsymbol{\beta}_k + \epsilon_{ki},
  }
where \eqn{\boldsymbol{X}^T\boldsymbol{\beta}} is a baseline mean function, and \eqn{\boldsymbol{X}} is high-dimensional.
}

\references{
Liu, X. (2024). High-dimensional test of relevant difference and its application to transfer learning. Manuscript.
}

\examples{
  data(simulData_test_gauss)
  y <- datahb$Y[1:5]
  dim(datahb$X)
  dim(datahb$Z)

  data(simulData_trans_gauss)
  y <- dataset[[1]]$Y
  dim(dataset[[1]]$X)

  dim(dataset[[2]]$X)
}

\keyword{datasets}
