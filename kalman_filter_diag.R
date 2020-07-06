kalman_filter_diag <- function(y, A, C, Q, R, init_x, init_V, model){
  
  # Kalman filter
  #
  # inputs:
  # y(:,t) - the observation at time t
  # A - the system matrix
  # C - the observation matrix 
  # Q - the system covariance 
  # R - the observation covariance
  # init_x - the initial state (column) vector 
  # init_V - the initial state covariance
  #
  # outputs (where X is the hidden state being estimated)
  # x(:,t) = E[X(:,t) | y(:,1:t)]
  # V(:,:,t) = Cov[X(:,t) | y(:,1:t)]
  # VV(:,:,t) = Cov[X(:,t), X(:,t-1) | y(:,1:t)] t >= 2
  
  os <- size(y)[1]
  T <- size(y)[2]
  ss <- size(A,1) 
  
  # set default params
  u <- c()
  B <- c()
  ndx <- c()
  
  x <- matrix(0L, nrow = ss, ncol = T)
  V <- (array(c(0),dim = c(ss,ss,T)))
  VV <- (array(c(0),dim = c(ss,ss,T)))
  loglik <- 0
  
  for (t in 1:T) {
    m <- model[t]
    if (t == 1) {
      prevx <- init_x
      prevV <- init_V
      initial <- 1
    } else {
      prevx <- x[ ,t-1,drop=FALSE]     
      prevV <- V[ , ,t-1]
      initial <- 0
    }
    
    result_kud <- kalman_update_diag(A[ , ,m], C[ , ,m], Q[ , ,m], R[ , ,m], 
                                     y[ ,t], prevx, prevV, initial)
    x[,t] <- result_kud$xnew
    V[,,t] <- result_kud$Vnew
    LL <- result_kud$loglik
    VV[,,t] <-  result_kud$VVnew
    loglik <- loglik + LL
  }
  
  return(list(x=x, V=V, VV=VV))
}
