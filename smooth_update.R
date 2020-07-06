smooth_update <- function(xsmooth_future, Vsmooth_future, xfilt, Vfilt, Vfilt_future, VVfilt_future, A, Q){
  
  # one step of the backwards RTS smoothing equations.
  #
  # inputs:
  # xsmooth_future = E[X_t+1|T]
  # Vsmooth_future = Cov[X_t+1|T]
  # xfilt = E[X_t|t]
  # Vfilt = Cov[X_t|t]
  # Vfilt_future = Cov[X_t+1|t+1]
  # VVfilt_future = Cov[X_t+1,X_t|t+1]
  # A = system matrix for time t+1
  # Q = system covariance for time t+1
  #
  # outputs:
  # xsmooth = E[X_t|T]
  # Vsmooth = Cov[X_t|T]
  # xpred = E[X(t+1) | t]
  
  xpred <- A%*%xfilt
  Vpred <- A%*%Vfilt%*%t(as.matrix(A))+Q
  J <- Vfilt%*%t(A)%*%pinv(Vpred) # smoother gain matrix
  xsmooth <- xfilt+J%*%(xsmooth_future-xpred)
  Vsmooth <- Vfilt+J%*%(Vsmooth_future-Vpred)%*%t(J)
  
  return(list(xsmooth=xsmooth, Vsmooth=Vsmooth))
}
