#' Define vector magnitude
#'
#' @param x1 a vector
#' @returnvector magnitude

Magnitude<-function(x1){
  Dimension<-length(x1)
  Entry<-rep(0,Dimension)
  for (i in 1:(Dimension)) {
    Entry[i]<-x1[i]^2
  }
  VectorMagnitude<-(sum(Entry))^(1/2)
  return(VectorMagnitude)
}
