\name{sasv8ch20sec20}
\alias{sasv8ch20sec20}
\docType{data}
\title{Sample data for canonical correlation analysis from the SAS manual
}
\description{This dataset includes six variables of fitness club data from the SAS manual.
}
\usage{data("sasv8ch20sec20")}

\details{
\describe{A 20x6 data matrix.
}
}

\source{
\url{https://v8doc.sas.com/sashtml/stat/chap20/sect20.htm}
}

\examples{
data(sasv8ch20sec20)

## Canonical Correlation Analysis
## X variables
X_vars <- c("Weight", "Waist", "Pulse")

## Y variables
Y_vars <- c("Chins", "Situps", "Jumps")

cancorr(X_vars=X_vars, Y_vars=Y_vars, data=sasv8ch20sec20)
}
