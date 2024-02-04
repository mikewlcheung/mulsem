\name{mulSEM-package}
\alias{mulSEM-package}
\alias{mulSEM}
\docType{package}
\title{
Some Multivariate Analyses using Structural Equation Modeling
}
\description{
A set of functions for some multivariate analyses utilizing a
structural equation modeling (SEM) approach through the 'OpenMx' package.
These analyses include canonical correlation analysis (CANCORR),
redundancy analysis (RDA), and multivariate principal component
regression (MPCR). It implements procedures discussed in
Gu and Cheung (2023) <doi:10.1111/bmsp.12301>, Gu, Yung, and Cheung
(2019) <doi:10.1080/00273171.2018.1512847>, and Gu et al. (2023)
<doi:10.1080/00273171.2022.2141675>.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>, Fei Gu <gu@vt.edu>,
  Yiu-Fai Yung <Yiu-Fai.Yung@sas.com>

Maintainer: Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\references{
Gu, F., & Cheung, M. W.-L. (2023). A Model-based approach to
multivariate principal component regression: Selection of principal
components and standard error estimates for unstandardized regression
coefficients. \emph{British Journal of Mathematical and Statistical
  Psychology}, \bold{76(3)}, 605-622. https://doi.org/10.1111/bmsp.12301

Gu, F., Yung, Y.-F., & Cheung, M. W.-L. (2019). Four covariance
structure models for canonical correlation analysis: A COSAN modeling
approach. \emph{Multivariate Behavioral Research}, \bold{54}(2),
192-223. https://doi.org/10.1080/00273171.2018.1512847
  
Gu, F., Yung, Y.-F., Cheung, M. W.-L. Joo, B.-K., & Nimon,
K. (2022). Statistical inference in redundancy analysis: A direct
covariance structure modeling approach. \emph{Multivariate Behavioral
  Research}, \bold{58(5)}, 877-893. https://doi.org/10.1080/00273171.2022.2141675 
}

\examples{
\donttest{
## Canonical Correlation Analysis
cancorr(X_vars=c("Weight", "Waist", "Pulse"),
        Y_vars=c("Chins", "Situps", "Jumps"),
        data=sas_ex1)

## Redundancy Analysis
rda(X_vars=c("x1", "x2", "x3", "x4"),
    Y_vars=c("y1", "y2", "y3"),
    data=sas_ex2)
	
## Multivariate Principal Component Regression	
mpcr(X_vars=c("AU", "CC", "CL", "CO", "DF", "FB", "GR", "MW"),
     Y_vars=c("IDE", "IEE", "IOCB", "IPR", "ITS"),
     pca="COR", pc_select=1,
     data=Nimon21)
}
}
