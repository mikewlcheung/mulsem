\name{rda}
\alias{rda}
\title{Redundancy Analysis
}
\description{It conducts a redundancy analysis using the OpenMx package.
  Missing data are handled with the full information maximum likelihood
  method when raw data are available. It provides standard errors on the
  standardized estimates.
}
\usage{
rda(X_vars, Y_vars, data=NULL, Cov, numObs, extraTries=50, ...)
}
\arguments{
  \item{X_vars}{A vector of characters of the X variables.}
  \item{Y_vars}{A vector of characters of the Y variables.}
  \item{data}{A data frame of raw data.}
  \item{Cov}{A covariance or correlation matrix if \code{data} is not
	available.}
  \item{numObs}{A sample size if \code{data} is not available.}
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
  Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2021). Statistical inference
  in redundancy analysis: A direct covariance structure modeling
  approach. \emph{Manuscript submited for publication}.
}  
\seealso{\code{\link[mulSEM]{Chittum19}}, \code{\link[mulSEM]{sasv8ch65sec26}}}

