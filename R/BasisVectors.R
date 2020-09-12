#' Generate Basis Vectors based on number of enzymes
#'
#' This function takes number of enzymes as input and
#' defines the orthonormal basis for the hyperplane perpendicular
#' to the vector <1,1,1...1> in Rn, where n is equal to the number
#' of enzymes

#'
#'@param NumberOfEnzymes the number of enzymes used in pathway
#'@return NormalizedVectors
#'@export

BasisVectors<-function(NumberOfEnzymes){
  #The first piece of code generates the vectors that will be fed into the Gram-Schmidt algorithm.
  MatrixFilling<-rep(0,NumberOfEnzymes*(NumberOfEnzymes-1))
  Vectors<-matrix(MatrixFilling,nrow=NumberOfEnzymes)
  for (i in 1:(NumberOfEnzymes-1)) {
    Vectors[i+1,i]<-(-i)
    Vectors[1:i,i]<-(1)
  }
  NormalizedVectors<-matrix(MatrixFilling,nrow=NumberOfEnzymes)
  for (j in 1:(NumberOfEnzymes-1)) {
    NormalizedVectors[1:NumberOfEnzymes,j]<-Vectors[1:NumberOfEnzymes,j]/Magnitude(Vectors[1:NumberOfEnzymes,j])
  }
  return(NormalizedVectors)
}

