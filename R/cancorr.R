cancorr <- function(X_vars, Y_vars, data=NULL, Cov, numObs,
                    model=c("CORR-W", "CORR-L", "COV-W", "COV-L"),
                    extraTries=50, ...) {

    model <- match.arg(model)
    switch(model,
           "CORR-W" = {analysis.cor=TRUE; type="W"}, 
           "CORR-L" = {analysis.cor=TRUE; type="L"},
           "COV-W" =  {analysis.cor=FALSE; type="W"},
           "COV-L" =  {analysis.cor=FALSE; type="L"})
           
    if (is.null(data)) {
        ## Cov or Cor as inputs
        mxdata <- mxData(observed=Cov, type="cov", numObs=numObs)
        
        mx.model <- mxModel("CanCorr", mxdata, mxFitFunctionML(),
                            mxExpectationNormal(covariance="expCov",
                                                dimnames=c(X_vars, Y_vars)))        
    } else {
        ## Raw data as inputs
        ## Sample Cov as starting values if raw data are given
        Cov <- cov(data, use="pairwise.complete.obs")

        mxdata <- mxData(observed=data, type="raw")
        
        mx.model <- mxModel("CanCorr", mxdata, mxFitFunctionML(),
                            mxExpectationNormal(covariance="expCov",
                                                means="expMean",
                                                dimnames=c(X_vars, Y_vars)))    
    }
  
    ## No. of X variables
    p <- length(X_vars)

    ## No. of Y variables
    q <- length(Y_vars)

    ## Make sure that p is NOT longer than q
    d  <- q-p

    if (d<0) stop("Length of 'Y_vars' cannot be shorter than length of 'X_vars'. Please swap the variables between them.\n")

    if (analysis.cor) {
        R <- cov2cor(Cov)
    } else {
        R <- Cov
    }
        
    ## Calculate the starting values
    Rxx <- R[X_vars, X_vars]
    Ryy <- R[Y_vars, Y_vars]
    Rxy <- R[X_vars, Y_vars]

    ## ei <- fda::geigen(Rxy, Rxx, Ryy)
    ## names(ei) <- c("cor", "xcoef", "ycoef")

    ## https://stats.stackexchange.com/questions/363425/raising-a-variance-covariance-matrix-to-a-negative-half-power
    ## https://stat.ethz.ch/pipermail/r-help/2008-April/160647.html
    "%^%" <- function(x, n) with(eigen(x), vectors %*% (values^n * t(vectors)))
    Rx.5 <- Rxx %^% (-0.5)
    Ry.5 <- Ryy %^% (-0.5)
    SVD <- svd( Rx.5 %*% Rxy %*% Ry.5)

    ## Checking
    ## all.equal( Rx.5 %*% Rxy %*% Ry.5,
    ##            SVD$u %*% diag(SVD$d) %*% t(SVD$v) )

    ## Starting values for A1
    if (type=="W") {
        A1.st <- Rx.5 %*% SVD$u
    } else {
        A1.st <- solve(Rx.5 %*% SVD$u)
    }

    ## Starting values for A2
    if (type=="W") {
        A2.st <- Ry.5 %*% SVD$v
    } else {
        if (d==0) {
            A2.st <- solve(Ry.5 %*% SVD$v)
        } else {
            A2.st <- matrix(rnorm(q*p, mean=0.3, sd=1), nrow=q, ncol=p)
        }
    }
        
    A1 <- mxMatrix("Full", nrow=p, ncol=p, free=TRUE,
                  values=A1.st,
                  labels=outer(1:p, 1:p, function(x, y) paste0("A1_", x, "_", y)),
                  name="A1")

    A2 <- mxMatrix("Full", nrow=q, ncol=p, free=TRUE,
                  values=A2.st,
                  labels=outer(1:q, 1:p, function(x, y) paste0("A2_", x, "_", y)),
                  name="A2")
    
    ## pxp matrix of zero
    p0p <- mxMatrix("Zero", nrow=p, ncol=p, name="p0p")
  
    ## Identity matrices in dimensions of p and d
    Id <- mxMatrix("Iden", nrow=d, ncol=d, name="Id") 
    Ip <- mxMatrix("Iden", nrow=p, ncol=p, name="Ip")
    
    if (d==0) {
        F <- mxAlgebra(rbind(cbind(t(A1), p0p),
                             cbind(p0p, t(A2))), name="F")
        mx.model <- mxModel(mx.model, F, A1, A2, p0p)       
    } else {

        if (d==1) {        
            A3 <- mxMatrix("Full", nrow=q, ncol=1, free=TRUE,
                           values=rnorm(q, mean=0.3, sd=1),
                           labels=outer(1:q, 1, function(x, y) paste0("A3_", x, "_", y)),
                           name="A3")
        } else {
            ## D_s: A dxd matrix of T/F elements.
            ## The last d_s=d*(d-1)/2 elements are arbitrarily fixed at 0
            d_s <- d*(d-1)/2
            D_s <- matrix(rep(c(TRUE, FALSE), times=c(q*d-d_s, d_s)), nrow=q, ncol=d)
            startvalues <- matrix(c(rnorm(q*d-d_s, mean=0.3, sd=1), rep(0, d_s)),
                                  nrow=q, ncol=d)
             
            A3 <- mxMatrix("Full", nrow=q, ncol=d, free=D_s, values=startvalues,
                           outer(1:q, 1:d, function(x, y) paste0("A3_", x, "_", y)),
                           name="A3")   
        }
        p0q <- mxMatrix("Zero", nrow=p, ncol=q, name="p0q")  
        d0p <- mxMatrix("Zero", nrow=d, ncol=p, name="d0p")
        F <- mxAlgebra(rbind(cbind(t(A1), p0q),
                             cbind(p0p, t(A2)),
                             cbind(d0p, t(A3))),
                       name="F")
        mx.model <- mxModel(mx.model, F, A1, A2, A3, p0p, p0q, d0p) 
    }

    ## Create Lambda matrix
    Lambda <- mxMatrix("Diag", nrow=p, ncol=p, free=TRUE, labels=paste0("l", 1:p),
                       values=SVD$d,
                       lbound=1e-10,
                       name="Lambda")

    ## Create a repeated contrast
    cont  <- cbind(diag(p-1), matrix(0, ncol=1, nrow=(p-1))) +
        cbind(matrix(0, ncol=1, nrow=(p-1)), (-1)*diag(p-1))
    Cont <- mxMatrix("Full", nrow=p-1, ncol=p, free=FALSE,
                     values=cont, name="Cont")

    ## (p-1)x1 zeros
    p_1_0_1  <- mxMatrix("Zero", nrow=p-1, ncol=1, name="p_1_0_1")

    ## Make sure that lambdas are in descending order
    Constraint_lambda <- mxConstraint(Cont %*% diag2vec(Lambda) > p_1_0_1,
                                      name="Constraint_lambda")
    mx.model <- mxModel(mx.model, Lambda, Cont, p_1_0_1, Constraint_lambda)
    
    if (d==0) {
        P <- mxAlgebra(rbind(cbind(Ip, Lambda),
                             cbind(Lambda, Ip)),
                       name="P")
        mx.model <- mxModel(mx.model, P, Ip)
    } else {
        d0p <- mxMatrix("Zero", nrow=d, ncol=p, name="d0p")
        P <- mxAlgebra(rbind(cbind(Ip, Lambda, t(d0p)),
                             cbind(Lambda, Ip, t(d0p)),
                             cbind(d0p, d0p,Id)),
                       name="P")
        mx.model <- mxModel(mx.model, P, Ip, Id, d0p)
    } 

    ## Analysis of correlation matrix
    if (analysis.cor) {
        ## Matrix to fit SDs of the variables
        SD <- mxMatrix("Diag", nrow=(p+q), ncol=(p+q), free=TRUE, 
                   values=sqrt(diag(Cov[c(X_vars, Y_vars), c(X_vars, Y_vars)])), 
                   labels=c(paste0("sdx", 1:p), paste0("sdy", 1:q)),
                   name="SD")

        ## (p+q) column vector of ones
        p_qOne1 <- mxMatrix("Unit", nrow=p+q, ncol=1, name="p_qOne1") 
        
        if (type=="L") {
            expCov <- mxAlgebra( SD %&% (F %&% P), name="expCov" )
            constraint <- mxConstraint( diag2vec(F %&% P) == p_qOne1, 
                                       name = 'constraint' )
        } else {
            ## type=="W"
            expCov <- mxAlgebra( SD %&% (solve(F) %&% P), name="expCov" )
            constraint <- mxConstraint( diag2vec(solve(F) %&% P) == p_qOne1, 
                                       name = 'constraint' )
        }
        mx.model <- mxModel(mx.model, expCov, SD, p_qOne1, constraint)
    ## Analysis of covariance matrix
    } else {
        if (type=="L") {
            expCov <- mxAlgebra( F %&% P, name="expCov" )
        } else {
            expCov <- mxAlgebra( solve(F) %&% P, name="expCov" )
        }
        mx.model <- mxModel(mx.model, expCov)
    }
    
 
    ## Add expMean structure if raw data are used as inputs
    if (!is.null(data)) {  
        expMean <- mxMatrix("Full", nrow=1, ncol=p+q, free=TRUE, 
                            values=apply(data[, c(X_vars, Y_vars)], 2, mean, na.rm=TRUE),
                            labels = c(paste0("mux_", seq_len(p)), paste0("muy_", seq_len(q))),
                            name="expMean")        
        mx.model <- mxModel(mx.model, expMean)   
    }

    if (extraTries==0) {
        mx.fit <- suppressMessages(mxRun(mx.model, ...))
    } else {
        mx.fit <- mxTryHard(mx.model, extraTries = extraTries, ...) 
        ## Run it one more time to minimize error
        mx.fit <- suppressMessages(mxRun(mx.fit))
    }
    
    #### Constraints for checking
    ## Diagonals should be 1: Equation (14)
    if (analysis.cor) {
        if (type=="L") {
            Constraint <- t(mxEval(diag2vec(F %&% P), mx.fit))
        } else {
            Constraint <- t(mxEval(diag2vec(solve(F) %&% P), mx.fit))
        }
        colnames(Constraint) <- c(X_vars, Y_vars)
    } else {
        Constraint <- NULL
    }
        
    ## A1 matrix for X_vars
    A1_est <- mxEval(A1, mx.fit)
    #dimnames(A1_est) <- list(X_vars, X_vars)

    ## SE of A1 matrix
    A1_SE <- suppressMessages(mxSE(A1, mx.fit))
    #dimnames(A1_SE) <- list(X_vars, X_vars)

    ## A2 matrix for Y_vars
    A2_est <- mxEval(A2, mx.fit)
    #dimnames(A2_est) <- list(Y_vars, Y_vars)

    ## SE of A2 matrix
    A2_SE <- suppressMessages(mxSE(A2, mx.fit))
    #dimnames(A2_SE) <- list(Y_vars, Y_vars)
    
    ## A3 matrix
    if (d>0) {
        A3_est <- mxEval(A3, mx.fit)
        A3_SE <- suppressMessages(mxSE(A3, mx.fit))
    } else {
        A3_est <- NULL
        A3_SE <- NULL
    }

    ## Lambdas
    Lambdas_est <- mxEval(diag(Lambda), mx.fit)
    Lambdas_SE <- suppressMessages(matrix(diag(mxSE(Lambda, mx.fit)), ncol=1))
   
    out <- list(mx.fit=mx.fit, mx.model=mx.model, model=model, Constraint=Constraint,
                A1_est=A1_est, A2_est=A2_est, A3_est=A3_est,
                A1_SE=A1_SE, A2_SE=A2_SE, A3_SE=A3_SE, Lambdas_est=Lambdas_est,
                Lambdas_SE=Lambdas_SE)  
    class(out) <- "CanCorr"
    out
}

print.CanCorr <- function(x, digits=4, ...) {
    if (!is.element("CanCorr", class(x)))
        stop("\"x\" must be an object of class \"CanCorr\".")

    cat("Model is", x$model, ".\n")
    
    if (x$model %in% c("CORR-W", "CORR-L")) {
        cat("Please check the constraint before interpreting the results.\n")
        cat("Constraint: The followings should be close to 1: ", x$Constraint, ".\n")
    }
      
    cat("\nA1 matrix (X_vars):\n")
    .mprint(x$A1_est, digits=digits)
    cat("\nA1 matrix (SE):\n")
    .mprint(x$A1_SE, digits=digits)
  
    cat("\nA2 matrix (Y_vas):\n")
    .mprint(x$A2_est, digits=digits)
    cat("\nA2 matrix (SE):\n")
    .mprint(x$A2_SE, digits=digits)

    if (!is.null(x$A3_est)) {
         cat("\nA3 matrix:\n")
         .mprint(x$A3_est, digits=digits)
         cat("\nA3 matrix (SE):\n")
         .mprint(x$A3_SE, digits=digits)  
    }
     
    cat("\nLambdas:\n")
    .mprint(x$Lambdas_est, digits=digits)
    
    cat("\nLambdas (SE):\n")
    .mprint(x$Lambdas_SE, digits=digits)   
}
