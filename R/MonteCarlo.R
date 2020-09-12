#' MonteCarlo
#'
#' Function that generates the starting points to use for gradient descent,
#' only allowing values that perform better than an initial guess.
#' The MCGD algorithm is a combination of the MC and the GradientDescent functions.
#'
#' @param NumberOfEnzymes
#' @param TotalEnzyme
#' @param NumberOfRuns
#' @param ThresholdGuess
#' @param GDPrecision
#' @param GDNumberOfSteps
#' @return starting points for gradient descent algorithm


MC<-function(NumberOfEnzymes,TotalEnzyme,NumberOfRuns,ThresholdGuess,GDPrecision,GDNumberOfSteps){
  TotalEnzyme<-TotalEnzyme
  Threshold<-OFFunction(ThresholdGuess)
  AcceptableStarts<-0
  StartsNeeded<-NumberOfRuns
  StartMatrixFilling<-rep(0,(NumberOfEnzymes+1)*(NumberOfRuns))
  StartMatrix<-matrix(StartMatrixFilling,nrow=(NumberOfEnzymes+1))
  repeat{
    Total<-TotalEnzyme
    Random<-rep(0,NumberOfEnzymes)
    i<-0
    repeat{
      i<-i+1
      Random[i]<-runif(1,min=0,max=1)*Total
      Total<-Total-Random[i]
      if(i==(NumberOfEnzymes-1)){break}
    }
    Random[NumberOfEnzymes]<-Total
    Integers<-c(1:NumberOfEnzymes)
    OrderGenerator<-sample(Integers,NumberOfEnzymes,replace=FALSE)
    TestValueVector<-rep(0,NumberOfEnzymes)
    for (j in 1:NumberOfEnzymes) {
      TestValueVector[j]<-Random[OrderGenerator[j]]
    }
    TestValue<-OFFunction(TestValueVector)
    if(TestValue<Threshold){
      AcceptableStarts<-AcceptableStarts+1
      for (k in 1:NumberOfEnzymes) {
        StartMatrix[k,AcceptableStarts]<-TestValueVector[k]
      }
      StartMatrix[NumberOfEnzymes+1,AcceptableStarts]<-TestValue
    }
    if(AcceptableStarts==StartsNeeded){break}
  }
  return(StartMatrix)
}

