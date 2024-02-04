[![R build status](https://github.com/mikewlcheung/mulsem/workflows/R-CMD-check/badge.svg)](https://github.com/mikewlcheung/mulsem/actions)
[![cran version](http://www.r-pkg.org/badges/version/mulSEM)](https://cran.r-project.org/package=mulSEM)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/mulSEM)](https://cranlogs.r-pkg.org/badges/mulSEM)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/mulSEM)](https://cranlogs.r-pkg.org/badges/grand-total/mulSEM)
[![Rdoc](http://www.rdocumentation.org/badges/version/mulSEM)](https://www.rdocumentation.org/packages/mulSEM)

The `mulSEM` package includes some multivariate analyses utilizing a structural equation modeling (SEM) approach through the 'OpenMx' package. These analyses include canonical correlation analysis (CANCORR), redundancy analysis (RDA), and multivariate principal component regression (MPCR).

You may install it from CRAN by:

```
install.packages("mulSEM")
```

The developmental version can be installed from GitHub by:
```
## Install remotes package if it has not been installed yet
# install.packages("remotes")

remotes::install_github("mikewlcheung/mulsem")

library(mulSEM)

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
```

