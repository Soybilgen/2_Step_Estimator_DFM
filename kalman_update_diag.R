kalman_update_diag <- function(A, C, Q, R, y, x, V, initial){
  
  # do a one step update of the Kalman filter
  #
  # inputs:
  # A - the system matrix
  # C - the observation matrix 
  # Q - the system covariance 
  # R - the observation covariance
  # y(:)   - the observation at time t
  # x(:) - E[X | y(:, 1:t-1)] prior mean
  # V(:,:) - Cov[X | y(:, 1:t-1)] prior covariance
  #
  # outputs (where X is the hidden state being estimated)
  #  xnew(:) =   E[ X | y(:, 1:t) ] 
  #  Vnew(:,:) = Var[ X(t) | y(:, 1:t) ]
  #  VVnew(:,:) = Cov[ X(t), X(t-1) | y(:, 1:t) ]
  #  loglik = log P(y(:,t) | y(:,1:t-1)) log-likelihood of innovations
  
  # set default params
  u <- c()
  B <- c()
  
  if (initial==1) {
    if (isempty(u)) {
      xpred <- x
    } else { 
      xpred <- x+(B%*%u)
    }
    Vpred <- V
  } else {
    if (isempty(u)){
      xpred <- A%*%x
    } else{
      xpred <- (A%*%x)+(B%*%u)
    }
    Vpred <- (A%*%V%*%t(A))+Q
  }
  
  e <- y-C%*%xpred  # error (innovation)
  n <- max(size(e))
  ss <- max(size(A))
  d <- size(e,1)
  S <- (C %*% Vpred %*% t(C)) + R
  GG <- t(C)%*%diag(1/diag(R))%*%C
  Sinv <- diag(1/diag(R))-diag(1/diag(R))%*%C%*%pinv(diag(ss)+(Vpred%*%GG))%*%Vpred%*%t(C)%*%diag(1/diag(R))
  detS <- prod(diag(R))%*%det(diag(ss)+Vpred%*%GG)
  denom <- (2*pi)^(d/2)*sqrt(abs(detS))
  mahal <- apply(t(e)%*%Sinv%*%e,1, sum)
  loglik = -0.5%*%mahal - log(denom)
  K <- Vpred%*% t(C)%*%Sinv; # Kalman gain matrix
  
  xnew <- xpred+K%*%e               
  Vnew <- (diag(ss)-K%*%C)%*%Vpred    
  VVnew <- (diag(ss)-K%*%C)%*%A%*%V
  
  return(list(xnew=xnew, Vnew=Vnew, loglik=loglik, VVnew=VVnew))
}
