\name{cancorr}
\alias{cancorr}
\title{Canonical Correlation Analysis
}
\description{It conducts a canonical correlation analysis using the
  OpenMx package. Missing data are handled with the full information maximum likelihood
  method when raw data are available. It provides standard errors on the estimates.
}
\usage{
cancorr(X_vars, Y_vars, data=NULL, Cov, numObs,
        model=c("CORR-W", "CORR-L", "COV-W", "COV-L"),
        extraTries=50, ...)
}
\arguments{
  \item{X_vars}{A vector of characters of the X variables.}
  \item{Y_vars}{A vector of characters of the Y variables.}
  \item{data}{A data frame of raw data.}
  \item{Cov}{A covariance or correlation matrix if \code{data} is not
	available.}
  \item{numObs}{A sample size if \code{data} is not available.}
  \item{model}{Four models defined in Gu, Yung, and Cheung
	(2019). \code{CORR} and \code{COV} refer to the analysis of
	correlation structure and covariance structure, respectively. }
  \item{extraTries}{This function calls \code{\link[OpenMx]{mxTryHard}}
	to obtain the parameter estimates and their standard
	errors. \code{extraTries} indicates the number of extra runs. If
	\code{extraTries=0}, \code{\link[OpenMx]{mxRun}} is called.}
   \item{\dots}{Additional arguments sent to either
	\code{\link[OpenMx]{mxTryHard}} or \code{\link[OpenMx]{mxRun}}.}
}
\note{\code{cancorr} expects that there are equal or more number of
  variables in \code{Y_vars}. If there are fewer variables in
  \code{Y_vars}, you may swap between \code{X_vars} and \code{Y_vars}.}
\value{A list of output. The mx fitted object is in the slot of \code{mx.fit}.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\references{
  Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2019). Four covariance
  structure models for canonical correlation analysis: A COSAN modeling
  approach. \emph{Multivariate Behavioral Research}, \bold{54(2)},
  192-223. https://doi.org/10.1080/00273171.2018.1512847
}  
\seealso{\code{\link[mulSEM]{Thorndike00}}, \code{\link[mulSEM]{sasv8ch20sec20}}}

