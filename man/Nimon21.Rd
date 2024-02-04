\name{Nimon21}
\alias{Nimon21}
\docType{data}
\title{Raw data used in Nimon, Joo, and Bontrager (2021)
}
\description{This dataset includes the raw data of 13 variables reported by Nimon, Joo, and Bontrager (2021).
}
\usage{data("Nimon21")}

\details{
A data frame of 13 variables.
}

\source{
  Nimon, K., Joo, B.-K. (Brian), & Bontrager, M. (2021). Work cognitions and work intentions: A canonical correlation study. \emph{Human Resource Development International}, \bold{24}(1), 65-91. https://doi.org/10.1080/13678868.2020.1775038
}
\references{
  Gu, F., & Cheung, M. W.-L. (2023). A Model-based approach to multivariate principal component regression: Selection of principal components and standard error estimates for unstandardized regression coefficients. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{76(3)}, 605-622. 
  https://doi.org/10.1111/bmsp.12301 
  Gu, F., Yung, Y.-F., Cheung, M. W.-L. Joo, B.-K., & Nimon, K. (2023). Statistical inference in redundancy analysis: A direct covariance structure modeling
approach. \emph{Multivariate Behavioral Research}, \bold{58(5)}, 877-893. https://doi.org/10.1080/00273171.2022.2141675
}
\examples{
\donttest{
data(Nimon21)

## Redundancy Analysis
rda(X_vars=c("AU", "CC", "CL", "CO", "DF", "FB", "GR", "MW"),
    Y_vars=c("IDE", "IEE", "IOCB", "IPR", "ITS"),
    data=Nimon21)

## Multivariate Principal Component Regression	
mpcr(X_vars=c("AU", "CC", "CL", "CO", "DF", "FB", "GR", "MW"),
     Y_vars=c("IDE", "IEE", "IOCB", "IPR", "ITS"),
     pca="COR", pc_select=1,
     data=Nimon21)
}
}
