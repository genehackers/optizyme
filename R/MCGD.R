#' Gradient Descent using Monte Carlo sampling of starting points
#'
#' Function that takes the same inputs as GradientDescent
#' as well as a total amount of enzyme to be used as
#' the constraint in the optimization process,
#' an initial guess for a starting point for gradient descent,
#' and the number of starting points to use.
#' The function incorporates Monte Carlo sampling of starting points,
#' and only accepts starting points that have lower values
#' than the value fo the initial guess value.
#' A default guess for the starting point can be
#' the enzymes expressed in a 1:1:1… ratio.
#'
#' @param NumberOfEnzymes integer number of enzymes used in pathway
#' @param TotalEnzyme total enzyme
#' @param NumberOfRuns integer number of runs
#' @param ThreshholdGuess a vector
#' @param GDPrecision a float number
#' @param GDNumberOfSteps number of steps taken
#' @return SolutionVector

#Before MCGD can be used, a function named OFFunction that defines the model and outputs what the algorithm will optimize with respect to must be run.
MCGD<-function(NumberOfEnzymes,TotalEnzyme,NumberOfRuns,ThresholdGuess,GDPrecision,GDNumberOfSteps,GDdtCutoff){
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
  #The chunk of code below will be the gradient descent portion of the algorithm that is conducted on the acceptable starts determined by the Monte Carlo above.
  SearchMatrix<-matrix(StartMatrixFilling,nrow=(NumberOfEnzymes+1))
  i<-0
  repeat{
    i<-i+1
    Object<-GradientDescent(NumberOfEnzymes,StartMatrix[1:NumberOfEnzymes,i],GDNumberOfSteps,GDPrecision,GDdtCutoff)
    for (a in 1:(NumberOfEnzymes+1)) {
      SearchMatrix[a,i]<-Object[a]
    }
    if(i==AcceptableStarts){break}
  }
  Minimum<-min(SearchMatrix[NumberOfEnzymes+1,1:StartsNeeded])
  for (b in 1:NumberOfRuns) {
    if(SearchMatrix[(NumberOfEnzymes+1),b]==Minimum){
      SolutionVector<-SearchMatrix[1:(NumberOfEnzymes+1),b]
    }
  }
  return(SolutionVector)
}
