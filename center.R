center <- function(x) {
  # centers each column of X
  T <- dim(x)[1]
  N <- dim(x)[2]
  xc <- x-(matrix(1, nrow = T, ncol = N)%*%diag(apply(x, 2, sum)/T))
  return(xc)
}