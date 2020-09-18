#' EnzymeModelPlot
#'
#' This function provides visualization of the time course for a cell-free system
#'
#' @param NumberOfEnzymes is a scalar value equal to the number of enzymes in the system
#' @param InitialState is a vector containing the initial concentrations of each fo the species involved in the system. The last entry in InitialState should be a zero, and won't represent any species in the system. This last entry being zero is used for indexing purposes by the function.
#' @param Time is a scalar value that represents how long the simulated time course should represent. The units for time will depend on the units of the inputs to the function. For example, if turnover rates of enzymes are expressed in inverse seconds, the units of time will be seconds.
#' @param dt is a scalar value that represents the time step used in the integration of the differential equations representing the enzyme time course. As dt decreases, the computational cost increases, but accuracy of the model also increases.
#' @param EnzymeNumberOfSubstrates is a vector of size n, where n is equal to the number of enzymes. For example, if the first enzyme is a single substrate enzyme, then the first entry in EnzymeNumberOfSubstrates should be 1. If the second enzyme in the pathway is a bisubstrate enzyme, then the second entry in EnzymeNumberOfSubstrates should be 2. It is important to note that these functions do not support more than two substrates per enzyme.
#' @param EnzymeSubstrates is a matrix of size 2xn, where n is equal to the number of enzymes. Its entries are which species in the InitialState vector are substrates to the nth enzyme in the pathway. For example, if the first enzyme has one substrate and it is compound 1 listed in InitialState, then the first column of the EnzymeSubstrates matrix is 1,0. There must be two entries into the matrix for each column, and if one of the slots isnt used then the value should be 0. The first entry in the column shoul d be filled first. If the second enzyme has two substrates and they are compounds 2 and 3, then the second column of EnzymeSubstrates should be 2,3.
#' @param EnzymeProducts is also a matrix of size 2xn, with the same structure and properties as EnzymeSubstrates. However, the entries should represent the products of the enzyme rather than the substrates.
#' @param EnzymeConcentrations is a vector of size n. The first entry in the vector should be the concentration of the first enzyme in the pathway. The second entry is the concentration of the second enzyme, and so on.
#' @param EnzymeTurnoverRates is a vector of size n. It should be filled in the same manner as EnzymeConcentrations.
#' @param EnzymeMichaelisConstants is a matrix of dimensions 2xn. The entries in the matrix are the corresponding Michaelis constants for the enzymatic catalysis of each of the substrates. For example, in our example for EnzymeSubstrates, we said that the first column of EnzymeSubstrates was 1,0. The corresopnding column of EnzymeMichaelisConstants should be Km,1, where Km is the Michaelis constant for the catalysis of compound 1 by enzyme 1 in this example. If a Michaelis Constant does not exist because an enzyme is single substrate, the second entry in the column should be a 1. However, any nonzero constant will work here.
#' @param EnzymeS1CompetitiveInhibitors is a matrix of size 10xn, where n is the number of enzymes in the system. Each column of the matrix contains the competitive inhibitors with the first substrate listed in the corresponding EnzymeSubstrates matrix. If there are unused slots, they should be filled in with the final entry of InitialState (Assuming that the final entry is the entry that doesn't represent a species and is just zero). For example, if InitialState is a vector of size 12, then unused slots in the columns of EnzymeS1CompetitiveInhibitors should have the value 12. If compounds 5,6,7 are competitive inhibitors of compound 1 for enzyme 1, then the first column of the matrix is <5,6,7,12,12,12,12,12,12,12>, where 12 is the terminal entry of InitialState.
#' @param EnzymeS1CompetitiveInhibitionConstants is a matrix of the same dimensions as EnzymeS1CompetitiveInhibitors, and contains the corresponding competitive inhibition constants for each of the competitive inhibitors with the first substrate binding site. The competitive inhibition constants associated with unused competitive inhibitors should be arbitrary nonzero constants, it is recommended to use zero.
#' @param EnzymeS2CompetitiveInhibitors is structurally the same as its S1 counterpart, but the inputs differ based on whether the compounds are competitive inhibitors with the first substrate or the second substrate listed in the EnzymeSubstrates matrix.
#' @param EnzymeS2CompetitiveInhibitionConstants  is structurally the same as its S2 counterpart, but the inputs differ based on whether the compounds are competitive inhibitors with the first substrate or the second substrate listed in the EnzymeSubstrates matrix.
#' @param EnzymeNonCompetitiveInhibitors is also a matrix of size 10xn and should be filled in the same manner as the competitive inhibition matrices.
#' @param EnzymeNonCompetitiveInhibitorConstants is also a matrix of size 10xn and should be filled in the same manner as the competitive inhibition matrices.
#' @param FinalProduct ask pat
#' @return a plot of the time course of the system

EnzymeModelPlot<-function(NumberOfEnzymes,InitialState,Time,dt,EnzymeNumberOfSubstrates,EnzymeSubstrates,EnzymeProducts,EnzymeConcentrations,EnzymeTurnoverRates,EnzymeMichaelisConstants,EnzymeS1CompetitiveInhibitors,EnzymeS1CompetitiveInhibitionConstants,EnzymeS2CompetitiveInhibitors,EnzymeS2CompetitiveInhibitionConstants,EnzymeNonCompetitiveInhibitors,EnzymeNonCompetitiveInhibitionConstants,FinalProduct){

  NumberOfSteps<-Time/dt
  NumberOfCompounds<-length(InitialState)
  TimeCourseMatrixFilling<-rep(0,NumberOfCompounds*NumberOfSteps)
  TimeCourseMatrix<-matrix(TimeCourseMatrixFilling,ncol=NumberOfSteps)
  TimeCourseMatrix[1:NumberOfCompounds,1]<-InitialState
  Counter<-0
  repeat{
    Counter<-Counter+1

    #The first thing we need to do is define the rate at which each of the enzymes work.
    EnzymeRateVector<-rep(0,NumberOfEnzymes)
    for (a in 1:NumberOfEnzymes) {

      if(EnzymeNumberOfSubstrates[a]==1){

        CompetitiveInhibitorConcentrations<-rep(0,10)
        for (i in 1:10) {
          CompetitiveInhibitorConcentrations[i]<-TimeCourseMatrix[EnzymeS1CompetitiveInhibitors[i,a],Counter]
        }

        #The piece of code below generates the vector of non-competitive inhibitor concentrations that can be used as input in the Rate1 function.
        NonCompetitiveInhibitorConcentrations<-rep(0,10)
        for (j in 1:10) {
          NonCompetitiveInhibitorConcentrations[j]<-TimeCourseMatrix[EnzymeNonCompetitiveInhibitors[j,a],Counter]
        }

        EnzymeRateVector[a]<-RateFunction1(EnzymeConcentrations[a],EnzymeTurnoverRates[a],TimeCourseMatrix[EnzymeSubstrates[1,a],Counter],EnzymeMichaelisConstants[1,a],CompetitiveInhibitorConcentrations,EnzymeS1CompetitiveInhibitionConstants[1:10,a],NonCompetitiveInhibitorConcentrations,EnzymeNonCompetitiveInhibitionConstants[1:10,a])
      }

      ###################################################################
      if(EnzymeNumberOfSubstrates[a]==2){
        E1S1CompetitiveInhibitorConcentrations<-rep(0,10)
        for (i in 1:10) {
          E1S1CompetitiveInhibitorConcentrations[i]<-TimeCourseMatrix[EnzymeS1CompetitiveInhibitors[i,a],Counter]
        }

        E1S2CompetitiveInhibitorConcentrations<-rep(0,10)
        for (j in 1:10) {
          E1S2CompetitiveInhibitorConcentrations[j]<-TimeCourseMatrix[EnzymeS2CompetitiveInhibitors[j,a],Counter]
        }

        E1NonCompetitiveInhibitorConcentrations<-rep(0,10)
        for (k in 1:10) {
          E1NonCompetitiveInhibitorConcentrations[k]<-TimeCourseMatrix[EnzymeNonCompetitiveInhibitors[k,a],Counter]
        }

        EnzymeRateVector[a]<-RateFunction2(EnzymeConcentrations[a],EnzymeTurnoverRates[a],TimeCourseMatrix[EnzymeSubstrates[1:2,a],Counter],EnzymeMichaelisConstants[1:2,a],E1S1CompetitiveInhibitorConcentrations,EnzymeS1CompetitiveInhibitionConstants[1:10,a],E1S2CompetitiveInhibitorConcentrations,EnzymeS2CompetitiveInhibitionConstants[1:10,a],E1NonCompetitiveInhibitorConcentrations,EnzymeNonCompetitiveInhibitionConstants[1:10,a])
      }
    }

    RateVector<-rep(0,length(InitialState))
    for (b in 1:NumberOfEnzymes) {

      if(EnzymeNumberOfSubstrates[b]==1){
        for (i in 1:length(InitialState)) {
          if(i==EnzymeSubstrates[1,b]){RateVector[i]<-RateVector[i]-EnzymeRateVector[b]}
          if(i==EnzymeProducts[1,b]){RateVector[i]<-RateVector[i]+EnzymeRateVector[b]}
        }}

      if(EnzymeNumberOfSubstrates[b]==2){
        for (i in 1:length(InitialState)) {
          if(i==EnzymeSubstrates[1,b]){RateVector[i]<-RateVector[i]-EnzymeRateVector[b]}
          if(i==EnzymeSubstrates[2,b]){RateVector[i]<-RateVector[i]-EnzymeRateVector[b]}
          if(i==EnzymeProducts[1,b]){RateVector[i]<-RateVector[i]+EnzymeRateVector[b]}
          if(i==EnzymeProducts[2,b]){RateVector[i]<-RateVector[i]+EnzymeRateVector[b]}
        }}

    }

    for (i in 1:length(InitialState)) {
      TimeCourseMatrix[i,Counter+1]<-TimeCourseMatrix[i,Counter]+RateVector[i]*dt
    }

    if(Counter==(NumberOfSteps-1)){break}
  }

  Integers<-(1:(NumberOfSteps))
  Time<-dt*Integers
  rainbow<-rainbow(length(TimeCourseMatrix[,1]))

  if(length(TimeCourseMatrix[,1])==3){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    legend("bottomright",legend=c("S1","S2","S3"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==4){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    legend("bottomright",legend=c("S1","S2","S3","S4"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==5){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    lines(Time,TimeCourseMatrix[5,],col=rainbow[5])
    legend("bottomright",legend=c("S1","S2","S3","S4","S5"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==6){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    lines(Time,TimeCourseMatrix[5,],col=rainbow[5])
    lines(Time,TimeCourseMatrix[6,],col=rainbow[6])
    legend("bottomright",legend=c("S1","S2","S3","S4","S5","S6"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==7){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    lines(Time,TimeCourseMatrix[5,],col=rainbow[5])
    lines(Time,TimeCourseMatrix[6,],col=rainbow[6])
    lines(Time,TimeCourseMatrix[7,],col=rainbow[7])
    legend("bottomright",legend=c("S1","S2","S3","S4","S5","S6","S7"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==8){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    lines(Time,TimeCourseMatrix[5,],col=rainbow[5])
    lines(Time,TimeCourseMatrix[6,],col=rainbow[6])
    lines(Time,TimeCourseMatrix[7,],col=rainbow[7])
    lines(Time,TimeCourseMatrix[8,],col=rainbow[8])
    legend("bottomright",legend=c("S1","S2","S3","S4","S5","S6","S7","S8"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==9){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    lines(Time,TimeCourseMatrix[5,],col=rainbow[5])
    lines(Time,TimeCourseMatrix[6,],col=rainbow[6])
    lines(Time,TimeCourseMatrix[7,],col=rainbow[7])
    lines(Time,TimeCourseMatrix[8,],col=rainbow[8])
    lines(Time,TimeCourseMatrix[9,],col=rainbow[9])
    legend("bottomright",legend=c("S1","S2","S3","S4","S5","S6","S7","S8","S9"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==10){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    lines(Time,TimeCourseMatrix[5,],col=rainbow[5])
    lines(Time,TimeCourseMatrix[6,],col=rainbow[6])
    lines(Time,TimeCourseMatrix[7,],col=rainbow[7])
    lines(Time,TimeCourseMatrix[8,],col=rainbow[8])
    lines(Time,TimeCourseMatrix[9,],col=rainbow[9])
    lines(Time,TimeCourseMatrix[10,],col=rainbow[10])
    legend("bottomright",legend=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10"),col=rainbow,title="Compounds",lty=.1,cex=.5)
  }

  if(length(TimeCourseMatrix[,1])==11){
    plot(Time,TimeCourseMatrix[1,],type="l",col=rainbow[1],xlab="Time",ylab="Concentration",main="Concentrations Over Time")
    lines(Time,TimeCourseMatrix[2,],col=rainbow[2])
    lines(Time,TimeCourseMatrix[3,],col=rainbow[3])
    lines(Time,TimeCourseMatrix[4,],col=rainbow[4])
    lines(Time,TimeCourseMatrix[5,],col=rainbow[5])
    lines(Time,TimeCourseMatrix[6,],col=rainbow[6])
    lines(Time,TimeCourseMatrix[7,],col=rainbow[7])
    lines(Time,TimeCourseMatrix[8,],col=rainbow[8])
    lines(Time,TimeCourseMatrix[9,],col=rainbow[9])
    lines(Time,TimeCourseMatrix[10,],col=rainbow[10])
    lines(Time,TimeCourseMatrix[11,],col=rainbow[11])
    legend("bottomright",legend=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11"),col=rainbow,title="Compounds",lty=1,cex=.5)
  }

}
