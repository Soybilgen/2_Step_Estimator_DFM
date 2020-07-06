FactorExtraction <- function(x, q, r, p) {
  
  # extract common factors from vector of time series possibly unbalanced
  # at the end of the sample
  
  # The model
  # x_t = C F_t + \xi_t
  # F_t = AF_{t-1} + B u_t
  # R = E(\xi_t \xi_t')
  # Q = BB'
  # u_t ~ WN(0,I_q)
  # initx = F_0
  # initV = E(F_0 F_0')
  # ss: std(x) 
  # MM: mean(x)
  
  # q: dynamic rank
  # r: static rank (r>=q)
  # p: ar order of the state vector (default p=1)
  
  # F: estimated factors
  # VF: estimation variance for the common factors
  
  if (r < q) {
    message('q has to be less or equal to r')
    stop()
  }
  
  if (p < 1) {
    message("p has to be greater or equal to 0")
    stop()
  }
  
  x <- as.matrix(x)
  
  library("pracma") # Numeric Analysis Package
  library("RSpectra") # Solvers for Large-Scale Eigenvalue and SVD Problems
  
  T <- dim(x)[1]  # dimension of the panel
  N <- dim(x)[2]
  
  # construct the balanced panel z from the original panel x
  das <- colSums(is.na(x))
  m <- max(das)
  selected_num_row <- T-m
  z <- x[1:selected_num_row, ]
  
  # standardize the panel
  ss <- apply(z, 2, sd)
  MM <- apply(z, 2, mean)
  s <- matrix(1, nrow = T, ncol = length(ss))%*%diag(ss)
  M <- matrix(1, nrow = T, ncol = length(ss))%*%diag(MM)
  x = (x - M)/s
  z <- x[1:selected_num_row, ]
  
  # Estimate the parameters using the balanced part of the panel
  result_ricsw <- ricSW(z, q, r, p)
  A <- result_ricsw$A
  C <- result_ricsw$C
  Q <- result_ricsw$Q
  R <- result_ricsw$R
  initx <- result_ricsw$initx
  initV <- result_ricsw$initV
  
  # the signal extraction in presence of missing data is performed by
  # using a time varying Kalman filter in which missing data are assigned an
  # extremely large variance of the noise in idiosyncratic component.
  AA <- (array(c(A),dim = c(dim(A),T)))
  QQ <- (array(c(Q),dim = c(dim(Q),T)))
  CC <- (array(c(C),dim = c(dim(C),T)))
  RR <- (array(c(R),dim = c(dim(R),T)))
  
  for (jt in 1:T){
    miss <- is.nan(x[jt,])
    Rtemp <- matrix(diag(R),ncol = 1)
    Rtemp[miss] <- 1e+32
    RR[,,jt] <- diag(c(Rtemp))
  }
  
  # missing data are assigned an arbitrary value...
  xx<- x
  xx[is.nan(x)] = 0 
  
  # run the kalman smoother on the time varying state space model
  model <- 1:T
  ksd_result <- kalman_smoother_diag(t(xx), AA, CC, QQ, RR, initx, initV, model)
  xsmooth <- ksd_result$xsmooth
  Vsmooth <- ksd_result$Vsmooth
  
  VF <- Vsmooth
  F <-  t(xsmooth)
  
  return(list(F=F, VF=VF))
}