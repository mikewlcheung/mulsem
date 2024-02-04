\name{sas_ex1}
\alias{sas_ex1}
\docType{data}
\title{Sample data for canonical correlation analysis from the SAS manual
}
\description{This dataset includes six variables of fitness club data from the SAS manual.
}
\usage{data("sas_ex1")}

\details{
\describe{A 20x6 data matrix.
}
}

\source{
\url{https://documentation.sas.com/doc/en/statcdc/14.2/statug/statug_cancorr_example01.htm}
}

\examples{
data(sas_ex1)

## Canonical Correlation Analysis
cancorr(X_vars=c("Weight", "Waist", "Pulse"),
        Y_vars=c("Chins", "Situps", "Jumps"),
        data=sas_ex1)
}
