#.mprint <- function(x, digits=4) print(round(x, digits), digits=digits)
.mprint <- function(x, digits=4) print(noquote(format(round(x, digits), digits=digits)))
