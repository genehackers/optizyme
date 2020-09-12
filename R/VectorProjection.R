#' Vector projection operation
#'
#' @param x1
#' @param x2
#' @return VectorProjection

Projection<-function(x1,x2){
  Term1<-DotProduct(x1,x2)
  Term2<-DotProduct(x2,x2)
  VectorProjection<-(Term1/Term2)*x2
  return(VectorProjection)
}
