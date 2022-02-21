\name{mpcr}
\alias{mpcr}
\title{Multivariate Principal Component Regression
}
\description{It conducts a multivariate principal component regression analysis using the OpenMx package. Missing data are handled with the full information maximum likelihood
  method when raw data are available. It provides standard errors on the
  estimates.
}
\usage{
mpcr(X_vars, Y_vars, data=NULL, Cov, Means=NULL, numObs, model=c("COV", "COR"), extraTries=50, ...)
}
\arguments{
  \item{X_vars}{A vector of characters of the X variables.}
  \item{Y_vars}{A vector of characters of the Y variables.}
  \item{data}{A data frame of raw data.}
  \item{Cov}{A covariance or correlation matrix if \code{data} is not
	available.}
  \item{Means}{An optional mean vector if \code{data} is not available.}
  \item{numObs}{A sample size if \code{data} is not available.}
  \item{model}{Whether the analysis is based unstandardized \code{COV}
	or standardized \code{COR}.}
  \item{extraTries}{This function calls \code{\link[OpenMx]{mxTryHard}}
	to obtain the parameter estimates and their standard
	errors. \code{extraTries} indicates the number of extra runs. If
	\code{extraTries=0}, \code{\link[OpenMx]{mxRun}} is called.}
   \item{\dots}{Additional arguments sent to either
	\code{\link[OpenMx]{mxTryHard}} or \code{\link[OpenMx]{mxRun}}.}
}

\value{A list of output. The mx fitted object is in the slot of \code{mx.fit}.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\references{
  Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2022). Standard error
  estimates for unstandardized regression coefficients in multivariate
  principal component regression. \emph{Manuscript in preparation}.
}  


