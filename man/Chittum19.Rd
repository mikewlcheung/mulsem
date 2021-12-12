\name{Chittum19}
\alias{Chittum19}
\docType{data}
\title{Correlation matrix of a model of motivation
}
\description{This dataset includes a correlation matrix of 12 variables
  (n=533) of a model of motivation reported by Chittum, Jones, and
  Carter (2019).
}
\usage{data("Chittum19")}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A 12x12 correlation matrix.}
  \item{n}{A sample size.}
}
}

\source{
Chittum, J. R., Jones, B. D., & Carter, D. M. (2019). A person-centered investigation of patterns in college students' perceptions of motivation in a course. \emph{Learning and Individual Differences}, \bold{69}, 94-107. https://doi.org/10.1016/j.lindif.2018.11.007
}
\references{
  Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2021). Statistical inference
  in redundancy analysis: A direct covariance structure modeling
  approach. \emph{Manuscript in preparation}.
}
\examples{
data(Chittum19)

## Redundancy Analysis
## X variables
X_vars <- c("Empowerment", "Usefulness", "Success", "Interest", "Caring")

## Y variables
Y_vars <- c("Final_Exam", "Learning", "Course_Rating", "Instr_Rating",
            "Effort", "Cog_Engage", "Cost")

rda(X_vars=X_vars, Y_vars=Y_vars, Cov=Chittum19$data, numObs=Chittum19$n)
}

