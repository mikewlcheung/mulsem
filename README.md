[![R build status](https://github.com/mikewlcheung/mulsem/workflows/R-CMD-check/badge.svg)](https://github.com/mikewlcheung/mulsem/actions)


The `mulSEM` package includes some multivariate analyses using a structural equation modeling (SEM) approach via the 'OpenMx' package. These analsyes include, e.g., canonical correlation analysis and redundancy analysis.

The developmental version can be installed from GitHub by:
```
## Install remotes package if it has not been installed yet
# install.packages("remotes")

remotes::install_github("mikewlcheung/mulsem")

library(mulSEM)

## Canonical Correlation Analysis
cancorr(X_vars=c("Weight", "Waist", "Pulse"),
        Y_vars=c("Chins", "Situps", "Jumps"),
        data=sasv8ch20sec20)

## Redundancy Analysis
rda(X_vars=c("x1", "x2", "x3", "x4"),
    Y_vars=c("y1", "y2", "y3"),
    data=sasv8ch65sec26)
```
