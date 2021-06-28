RA <- function(X_vars, Y_vars, data=NULL, Cov, numObs, extraTries=50, ...) {
  
  # oldw <- getOption("warn")
  # options(warn = -1)
  
  ## Sample Cov as starting values if raw data are given
  if (!is.null(data)) {
    Cov <- cov(data, use="pairwise.complete.obs")
  }
  
  ## No. of X variables
  p <- length(X_vars)

  ## No. of Y variables
  q <- length(Y_vars)
  
  ## Starting values for W
  ei <- eigen(solve(cov2cor(Cov[X_vars, X_vars]))) # eigen value decomposition
  start_W <- ei$vectors %*% diag(sqrt(ei$values)) 
  
  ## Starting values for Ly
  start_Ly <- cov2cor(Cov)[Y_vars, X_vars] %*% start_W
  
  ## Add a bit noise to the starting values
  start_W <- start_W + matrix(runif(p*p, min=-0.2, max=0.2), nrow=p, ncol=p)
  start_Ly <- start_Ly + matrix(runif(p*q, min=-0.2, max=.2), nrow=q, ncol=p)
  
  ## Standard deviations of X and Y
  SD <- mxMatrix("Diag", nrow=(p+q), ncol=(p+q), free=TRUE, 
                 values=sqrt(diag(Cov[c(X_vars, Y_vars), c(X_vars, Y_vars)])), 
                 labels=c(paste0("sdx", 1:p), paste0("sdy", 1:q)),
                 name="SD")
  
  ## qxp matrix of zero
  qZerop <- mxMatrix("Zero", nrow=q, ncol=p, name="qZerop")
  
  ## Identity matrices in dimensions of p and q
  Iq <- mxMatrix("Iden", nrow=q, ncol=q, name="Iq") 
  Ip <- mxMatrix("Iden", nrow=p, ncol=p, name="Ip") 

  ## W matrix with random starting values
  W <- mxMatrix("Full", nrow=p, ncol=p, free=TRUE,
                #values=rnorm(p*p, mean=0.3, sd=0.1),
                values=start_W,
                labels=outer(1:p, 1:p, function(x, y) paste0("W", x, "_", y)),
                name="W")

  ## Block diagonal matrix of W and I
  W_I <- mxAlgebra(rbind(cbind(solve(t(W)), t(qZerop)),
                         cbind(qZerop, Iq)), name="W_I")

  ## Ryy
  Ryy <- mxMatrix("Stand", nrow=q, ncol=q, free=TRUE,
                  values=vechs(cov2cor(Cov[Y_vars, Y_vars])),
                  labels=vechs(outer(1:q, 1:q, function(x, y) paste0("Ry", x, "_", y))),
                  name="Ryy")

  ## Ly
  if (p-q > 0) {
    ## Ly when p > q, the last (p-q) columns are zeros
    Ly <- mxMatrix("Full", nrow=q, ncol=p,
                   free=cbind(matrix(TRUE, nrow=q, ncol=q), 
                              matrix(FALSE, nrow=q, ncol=(p-q))),
                   #values=cbind(matrix(rnorm(q*q, mean=0.3, sd=0.1), nrow=q, ncol=q), 
                   #             matrix(0, nrow=q, ncol=(p-q))),
                   values=cbind(matrix(start_Ly[1:q, 1:q], nrow=q, ncol=q), 
                                matrix(0, nrow=q, ncol=(p-q))),
                   labels=outer(1:q, 1:p, function(x, y) paste0("Ly", x, "_", y)),
                   name="Ly")
  } else {
    ## Ly when p <= q
    Ly <- mxMatrix("Full", nrow=q, ncol=p,
                  free=TRUE,
                  #values=rnorm(p*q, mean=0.2, sd=0.1),
                  values=start_Ly,
                  labels=outer(1:q, 1:p, function(x, y) paste0("Ly", x, "_", y)),
                  name="Ly")
  }

  ## Matrix to store q, the no. of Y variables
  q_matrix <- mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=q, name="q")

  ## Lambdas: Equation (7)
  Lambdas <- mxAlgebra(diag2vec(t(Ly) %*% Ly)/q, name="Lambdas")

  ## Lx: Equation (12)
  Lx <- mxAlgebra(solve(t(W)), name="Lx")

  ## Block matrix of I and R
  I_R <- mxAlgebra(rbind(cbind(Ip, t(Ly)),
                         cbind(Ly, Ryy)), name="I_R")

  ## px1 column vector of ones
  pOne1 <- mxMatrix("Unit", nrow=p, ncol=1, name="pOne1")
  
  ## p*(p-1)/2x1 column vector of zeros 
  psZero1 <- mxMatrix("Zero", nrow=p*(p-1)/2, ncol=1, name="psZero1")

  ## Expected covariance matrix: Equation (10)
  expCov <- mxAlgebra(SD %&% (W_I %&% I_R), name="expCov")

  ## Equation (11)
  constraint1 <- mxConstraint( diag2vec(solve(t(W)) %*% solve(W)) == pOne1, 
                               name = 'constraint1' )

  ## Second constraint that off diagonals are zeros 
  constraint2 <- mxConstraint( vechs(t(Ly) %*% Ly) == psZero1, name = 'constraint2' )

  if (is.null(data)) {
    ## Cov or Cor as inputs
    mxdata <- mxData(observed=Cov[c(X_vars, Y_vars), c(X_vars, Y_vars)], 
                     type="cov", numObs=numObs)
    
    ## Expected covariance, which is required with raw data, in the fit function
    expFun <- mxExpectationNormal(covariance="expCov", dimnames=c(X_vars, Y_vars))
    
    ## Combine everything to a model
    model <- mxModel("RA", mxdata,
                     qZerop, SD, Iq, Ip, W, W_I, Ryy, Ly, I_R, pOne1, psZero1, Lambdas, 
                     Lx, q_matrix, expCov, expFun, mxFitFunctionML(), constraint1, constraint2)
  
  } else {
    ## Raw data as inputs
    mxdata <- mxData(observed=data[, c(X_vars, Y_vars)], type="raw" )
    
    expMean <- mxMatrix("Full", nrow=1, ncol=p+q, free=TRUE, 
                        values=apply(data, 2, mean, na.rm=TRUE),
                        labels = c(paste0("mux_", seq_len(p)), paste0("muy_", seq_len(q))),
                        name="expMean")
    
    ## Expected covariance in the fit function
    expFun <- mxExpectationNormal(covariance="expCov", means="expMean", 
                                  dimnames=c(X_vars, Y_vars))    
    
    ## Combine everything to a model
    model <- mxModel("RA", mxdata,
                     qZerop, SD, Iq, Ip, W, W_I, Ryy, Ly, I_R, pOne1, psZero1, 
                     Lambdas, expMean, Lx, q_matrix, expCov, expFun, 
                     mxFitFunctionML(), constraint1, constraint2)    
  }
  
  if (extraTries>0) {
    mx.fit <- mxTryHard(model, extraTries = extraTries, ...) 
  } else {
    mx.fit <- mxRun(model)
  }
  
  #### Constraints for checking
  ## Diagonals should be 1: Equation (14)
  Constraint1 <- mxEval(t(diag(Lx %*% t(Lx))), mx.fit)
  colnames(Constraint1) <- X_vars

  ## Off-diagonals should be 0
  Constraint2 <- mxEval(vechs(t(Ly) %*% Ly), mx.fit)
  Constraint2 <- c(min(Constraint2), max(Constraint2))
  names(Constraint2) <- c("Min", "Max")

  ## W matrix
  W_est <- mxEval(W, mx.fit)
  dimnames(W_est) <- list(X_vars, X_vars)

  ## Ly matrix
  Ly_est <- mxEval(Ly, mx.fit)
  dimnames(Ly_est) <- list(Y_vars, X_vars)

  ## Lx matrix
  Lx_est <- mxEval(Lx, mx.fit)
  dimnames(Lx_est) <- list(X_vars, X_vars)

  ## Lambdas
  Lambdas <- mxEval(Lambdas, mx.fit)

  ## SE of W matrix
  W_SE <- mxSE(W, mx.fit)
  dimnames(W_SE) <- list(X_vars, X_vars)
  
  ## SE of Ly matrix
  Ly_SE <- mxSE(Ly, mx.fit)
  dimnames(Ly_SE) <- list(Y_vars, X_vars)    
  
  ## SE of Lx matrix
  Lx_SE <- mxSE(Lx, mx.fit)
  dimnames(Lx_SE) <- list(X_vars, X_vars)   
    
  out <- list(mx.fit=mx.fit, Constraint1=Constraint1, Constraint2=Constraint2,
              W_est=W_est, Ly_est=Ly_est, Lx_est=Lx_est, Lambdas=Lambdas,
              W_SE=W_SE, Ly_SE=Ly_SE, Lx_SE=Lx_SE)  
  class(out) <- "RA"
  
  # options(warn = oldw)
  out
}

print.RA <- function(x, ...) {
  if (!is.element("RA", class(x)))
    stop("\"x\" must be an object of class \"RA\".")

  cat("Please check the constraints before interpreting the results.\n")
  cat("Constraint 1: The followings should be close to 1. ", x$Constraint1, "\n")
  cat("Constraint 2: The followings (min and max) should be close to 0. ", x$Constraint2, "\n")
  
  cat("\nW matrix:\n")
  print(x$W_est)
  cat("\nW matrix (SE):\n")
  print(x$W_SE)
  
  cat("\nLy matrix:\n")
  print(x$Ly_est)
  cat("\nLy matrix (SE):\n")
  print(x$Ly_SE)
  
  cat("\nLx matrix:\n")
  print(x$Lx_est)
  cat("\nLx matrix (SE):\n")
  print(x$Lx_SE)
  
  cat("\nLambdas:\n")
  print(x$Lambdas)
}
