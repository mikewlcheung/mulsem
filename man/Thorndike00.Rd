\name{Thorndike00}
\alias{Thorndike00}
\docType{data}
\title{Correlation matrix of a model of disgust
}
\description{This dataset includes a correlation matrix of 13 variables
  (n=679) between five subscales (y1 to y5) of the Disguest Emotion
  Scale and eight subscales (x1 to x8) of the Disgust Scale reported by
  Thorndike (2000, p. 238).
}
\usage{data("Thorndike00")}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A 13x13 correlation matrix.}
  \item{n}{A sample size.}
}
}

\source{
Thorndike, R. M. (2000). Canonical correlation analysis. In
 H. E. A. Tinsley & S. D. Brown (Eds.), \emph{Handbook of applied multivariate statistics and mathematical modeling}
(pp. 237-263). San Diego, CA: Academic Press.
}
\references{
Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2019). Four covariance structure models for canonical correlation analysis: A COSAN modeling approach. \emph{Multivariate Behavioral Research}, \bold{54(2)}, 192-223. https://doi.org/10.1080/00273171.2018.1512847
}
\examples{
\donttest{
data(Thorndike00)

## Canonical Correlation Analysis
## Note. We swap the X_vars and Y_vars because cancorr() expects that
## X_vars cannot have more variables than Y_vars.

cancorr(X_vars=c("y1", "y2", "y3", "y4", "y5"),
        Y_vars=c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"),
        Cov=Thorndike00$data, numObs=Thorndike00$n)
}
}
