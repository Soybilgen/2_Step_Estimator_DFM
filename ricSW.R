ricSW <- function(x, q, r, p) {
  
  # this function computes the parameters of the factor models
  
  # initialize
  Mx <- apply(x, 2, mean)
  Wx <- diag(apply(x, 2, sd)) 
  x <- center(x)%*%inv(Wx) 
  T <- dim(x)[1]
  N <- dim(x)[2]
  nlag <- p-1 # p=1, so nlag = 0.
  
  # define some preliminary quantity that are necessary to write the VAR in companion form
  A_temp <- matrix(0L, nrow = r, ncol = r*p)
  I <- diag(r*p)
  LL <- dim(I)[1]
  if (p != 1) {
    A <- rbind(A_temp, I[1:(LL-r), ])
  } else {
    A <- rbind(t(A_temp), I[0, ])
  }
  Q <- matrix(0L, nrow = r*p, ncol = r*p)
  Q[1:r, 1:r] <- diag(r)
  
  # computes eigenvalues and eigenvectors
  result_eigs <- eigs(cov(x), k=r, which = "LM")	
  # d is a rxr diagonal matrix with the largest eigenvalues on the diagonal 
  # v is a nxr matrix of the eigenvectors that corresponds to the eigenvalues
  d <- diag(length(result_eigs$values))*result_eigs$values
  v <- result_eigs$vectors
  v[,1]=v[,1]*-1
  F <- x%*%v  # PC estimates of the common factors
  R <- diag(diag(cov((x-x%*%v%*%t(v))))) # rstimate of the covariance matrix of the idiosincratic component
  
  # estimate the autoregressive model for the factors: run the var F(t) = A_1*F(t-1)+...+A_p*F(t-1) + e(t);
  z <- F
  Z <- c()
  for (kk in 1:p) {
    Z <- cbind(Z, z[(p-kk+1):(size(z)[1]-kk), ]) # stacked regressors (lagged SPC)
  }
  z <- z[(p+1):size(z)[1], ]
  A_temp <- (inv(t(Z)%*%Z)%*%t(Z))%*%z # OLS estimator of the VAR transition matrix
  A[1:r, 1:(r*p)] <- t(A_temp) 
  
  # compute Q
  e <- z-Z%*%A_temp # VAR residuals
  H <- cov(e) # VAR covariance matrix
  
  if (r==q){
    # the covariance matrix of the VAR residuals is of full rank
    Q[1:r, 1:r] = H
  } else {
    # the covariance matrix of the VAR residuals has reduced rank
    res_ed <- eigs(H, k=q, which = "LM") # eigenvalue decomposition
    P <- res_ed$vectors
    M <- res_ed$values
    M <- diag(length(M))*M
    P <- P%*%diag(sign(P[1, ]))
    Q[1:r, 1:r] = P%*%M%*%t(P)
  }
  
  # computes the initial conditions for the filter
  # the common factors are initialized by the PC estimates
  # initial variance is set equal to the unconditional variance of the common factors 
  z <- F
  Z <- c()
  for (kk in 0:nlag) {  
    Z <- cbind(Z, z[(nlag-kk+1):(size(z)[1]-kk), ]) # stacked regressors (lagged SPC)
  }
  initx <- t(t(Z[1,]))
  initV <- matrix((pinv(diag(size(kron(A,A),1))-kron(A, A))%*%(matrix(as.vector(Q), ncol = 1))), r*p, r*p)
  C <- cbind(v,matrix(0L, nrow = N, ncol = r*nlag)) 
  
  return(list(A=A, C=C, Q=Q, R=R, initx=initx, initV=initV, Mx=Mx, Wx=Wx))
}
