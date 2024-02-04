\name{sas_ex2}
\alias{sas_ex2}
\docType{data}
\title{Sample data for redundancy analysis from the SAS manual
}
\description{This dataset includes seven variables from the SAS manual.
}
\usage{data("sas_ex2")}

\details{
\describe{A 10x7 data matrix.
}
}

\source{
\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_transreg_details23.htm}
}

\examples{
data(sas_ex2)

## Redundancy Analysis
rda(X_vars=c("x1", "x2", "x3", "x4"),
    Y_vars=c("y1", "y2", "y3"),
    data=sas_ex2)
}
