\name{Lambert88}
\alias{Lambert88}
\docType{data}
\title{Correlation matrix of artificial data
}
\description{This dataset includes a correlation matrix of the
  artificial data 9 variables used in Table 1 of Lambert, Wildt, and
  Durand (1988).
}
\usage{data("Lambert88")}

\details{
A 9x9 correlation matrix.
}

\source{
Lambert, Z. V., Wildt, A. R., & Durand, R. M. (1988). Redundancy analysis: An alternative to canonical correlation and multivariate multiple regression in exploring interset associations. \emph{Psychological Bulletin}, \bold{104}(2), 282-289. https://doi.org/10.1037/0033-2909.104.2.282
}
\references{
Gu, F., Yung, Y.-F., Cheung, M. W.-L. Joo, B.-K., & Nimon, K. (2022). Statistical inference in redundancy analysis: A direct covariance structure modeling
approach. \emph{Multivariate Behavioral Research}. https://doi.org/10.1080/00273171.2022.2141675
}
\examples{
data(Lambert88)

## Redundancy Analysis
rda(X_vars=paste0("x", 1:5), Y_vars=paste0("y", 1:4), Cov=Lambert88, numObs=100)
}

