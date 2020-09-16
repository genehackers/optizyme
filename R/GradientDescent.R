#' Gradient Descent Script
#'
#' Function that takes the number of enzymes allowed to vary,
#' initial guess for the concentration of each of the
#' enzymes (inputted as a vector), number of steps to be performed,
#' and the multiple used to determine the step size of gradient descent.
#' Recommended that the multiple be below 1, however the algorithm will
#' adjust the step size after each iteration.
#'
#' @param NumberOfEnzymes the integer number of enzymes used in the pathway
#' @param InitialGuess a vector of inital guess values
#' @param NumberOfSteps an integer
#' @param GDdtMultiple a float number
#' @return a vector and the number of steps to get there

GradientDescent<-function(NumberOfEnzymes,InitialGuess,NumberOfSteps,GDdtMultiple,GDdtCutoff){
  NumberOfSteps<-NumberOfSteps
  GDMatrixFilling<-rep(0,(NumberOfEnzymes+1)*(NumberOfSteps))
  GDMatrix<-matrix(GDMatrixFilling,nrow=NumberOfEnzymes+1)
  GDMatrix[1:(NumberOfEnzymes),1]<-InitialGuess
  StartingConcentrations<-InitialGuess
  MinimumStarting<-min(StartingConcentrations)
  GDdt<-GDdtMultiple*MinimumStarting
  OrthonormalBasis<-BasisVectors(NumberOfEnzymes)
  Counter<-0
  repeat{
    Counter<-Counter+1
    if(Counter==NumberOfSteps){break}
    GDMatrix[NumberOfEnzymes+1,Counter]<-OFFunction(GDMatrix[1:NumberOfEnzymes,Counter])

    Adjacent<-rep(0,NumberOfEnzymes-1)
    for (i in 1:(NumberOfEnzymes-1)) {
      Adjacent[i]<-OFFunction(GDMatrix[1:NumberOfEnzymes,Counter]+OrthonormalBasis[1:NumberOfEnzymes,i]*GDdt)
    }


    Partial<-rep(0,NumberOfEnzymes-1)
    for (j in 1:(NumberOfEnzymes-1)) {
      Partial[j]<-(Adjacent[j]-GDMatrix[NumberOfEnzymes+1,Counter])/GDdt
    }


    GradientMagnitude<-Magnitude(Partial)
    if(GradientMagnitude==0){
      GDMatrix[1:(NumberOfEnzymes+1),NumberOfSteps]<-GDMatrix[1:(NumberOfEnzymes+1),Counter]
      break
    }


    GradientVector<-(Partial/GradientMagnitude)*GDdt
    TranslatedGradientFilling<-rep(0,NumberOfEnzymes*(NumberOfEnzymes-1))
    TranslatedGradient<-matrix(TranslatedGradientFilling,nrow=NumberOfEnzymes)
    for (k in 1:(NumberOfEnzymes-1)) {
      TranslatedGradient[1:NumberOfEnzymes,k]<-GradientVector[k]*OrthonormalBasis[1:NumberOfEnzymes,k]
    }


    ActualGradient<-rep(0,NumberOfEnzymes)
    for (a in 1:(NumberOfEnzymes)) {
      ActualGradient[a]<-sum(TranslatedGradient[a,1:(NumberOfEnzymes-1)])
    }



    GDMatrix[1:NumberOfEnzymes,Counter+1]<-GDMatrix[1:NumberOfEnzymes,Counter]-ActualGradient
    GDMatrix[(NumberOfEnzymes+1),Counter+1]<-OFFunction(GDMatrix[1:NumberOfEnzymes,Counter]-ActualGradient)

    if(GDMatrix[(NumberOfEnzymes+1),Counter+1]>GDMatrix[(NumberOfEnzymes+1),Counter]){
      GDdt<-GDdt/2
      GDMatrix[1:(NumberOfEnzymes+1),Counter+1]<-GDMatrix[1:(NumberOfEnzymes+1),Counter]
    }

    #The piece of code below is a kill switch when the accuracy of the algorithm becomes unnecessary.
    if(GDdt<GDdtCutoff*MinimumStarting){
      GDMatrix[1:(NumberOfEnzymes+1),NumberOfSteps]<-GDMatrix[1:(NumberOfEnzymes+1),Counter]
      break
    }

  }
  return(GDMatrix[1:(NumberOfEnzymes+1),NumberOfSteps])
}
StartingVector<-c(2,1.5,5.7,9.8,3.6)
GradientDescent(5,StartingVector,25,.1,.001)
