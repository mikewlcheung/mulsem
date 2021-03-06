\name{sasv8ch65sec26}
\alias{sasv8ch65sec26}
\docType{data}
\title{Sample data for redundancy analysis from the SAS manual
}
\description{This dataset includes seven variables from the SAS manual.
}
\usage{data("sasv8ch65sec26")}

\details{
\describe{A 10x7 data matrix.
}
}

\source{
\url{https://v8doc.sas.com/sashtml/stat/chap65/sect26.htm}
}

\examples{
data(sasv8ch65sec26)

## Redundancy Analysis
rda(X_vars=c("x1", "x2", "x3", "x4"),
    Y_vars=c("y1", "y2", "y3"),
    data=sasv8ch65sec26)
}
