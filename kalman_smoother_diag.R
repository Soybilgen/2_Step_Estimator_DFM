kalman_smoother_diag <- function(y, A, C, Q, R, init_x, init_V, model){
  
  # Kalman/RTS smoother.
  T <- size(y)[2]
  ss <- size(A,1)
  
  # set default params
  u <- c()
  B <- c()
  
  xsmooth <- matrix(0L, nrow = ss, ncol = T)
  Vsmooth <- (array(c(0),dim = c(ss,ss,T)))
  
  # forward pass
  kfd_result <- kalman_filter_diag(y, A, C, Q, R, init_x, init_V, model)
  xfilt  <- kfd_result$x
  Vfilt  <- kfd_result$V
  VVfilt <- kfd_result$VV
  
  # backward pass
  xsmooth[,T] = xfilt[,T]
  Vsmooth[,,T] = Vfilt[,,T]
  for (t in (T-1):1) {
    m <- model[t+1]
    result_s_update = smooth_update(xsmooth[ , t+1, drop=FALSE], as.matrix(Vsmooth[ , , t+1]), xfilt[ , t],
                                    as.matrix(Vfilt[ , , t]), as.matrix(Vfilt[ , , t+1]), as.matrix(VVfilt[ , , t+1]), 
                                    as.matrix(A[ , , m]), as.matrix(Q[ , , m]))
    xsmooth[,t] <- result_s_update$xsmooth
    Vsmooth[,,t] <- result_s_update$Vsmooth
  }
  
  return(list(xsmooth=xsmooth, Vsmooth=Vsmooth))
}
