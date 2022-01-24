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
cancorr(X_vars=c("Weight", "Waist", "Pulse"),
        Y_vars=c("Chins", "Situps", "Jumps"),
        data=sasv8ch20sec20)
}
