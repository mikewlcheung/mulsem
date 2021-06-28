\name{mulSEM-package}
\alias{mulSEM-package}
\alias{mulSEM}
\docType{package}
\title{
Some Multivariate Analyses using Structural Equation Modeling
}
\description{
  A collection of functions for some multivariate analyses using a
  structural equation modeling (SEM) approach via the 'OpenMx' package.
  These analsyes include, e.g., canonical correlation analysis and
  redundancy analysis, see Gu, Yung, and Cheung (2019)
  <doi:10.1080/00273171.2018.1512847>.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>, Fei Gu <gu@vt.edu>,
  Yiu-Fai Yung <Yiu-Fai.Yung@sas.com>

Maintainer: Mike Cheung <mikewlcheung@nus.edu.sg>
}
\references{
Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2019). Four covariance structure models for canonical correlation analysis: A COSAN modeling approach. \emph{Multivariate Behavioral Research}, \bold{54}(2), 192-223. https://doi.org/10.1080/00273171.2018.1512847
  
Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2021). Statistical inference
  in redundancy analysis: A direct covariance structure modeling
  approach. \emph{Manuscript in preparation}.

}

\examples{
## Seed for reproducibility
set.seed(12345)

## Redundancy Analysis
## X variables
X_vars <- c("Empowerment", "Usefulness", "Success", "Interest", "Caring")

## Y variables
Y_vars <- c("Final_Exam", "Learning", "Course_Rating", "Instr_Rating",
            "Effort", "Cog_Engage", "Cost")

RA(X_vars=X_vars, Y_vars=Y_vars, Cov=Chittum19$data, numObs=Chittum19$n)
}
