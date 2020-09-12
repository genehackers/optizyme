#'Defining dot product
#'
#'@param x1 a vector - explain
#'@param x2 a vector - explain
#'@return dot product to be used to generate basis vectors

DotProduct<-function(x1,x2){
  return(sum(x1*x2))
}
