#' example run of Monte Carlo Gradient Descent
#'
#' @return Solution Vector

example_optimization <- function(){
  OFFunction<-function(InitialGuess){
    NumberOfSteps<-6307
    dt<-.01
    VmXDH <- 123.697
    KmXDHXyl <- 0.19850960840785867
    KmXDHNAD <- 0.16123927045886804
    KmXDHXylLact <- 0.5357339888880199
    KmXDHNADH <- 0.02970905659171423
    kXLA <- 0.0071827237487217975
    VmXLA <- 839.404
    KmXLAXylLact <- 0.4452685845967363
    VmXAD <- 43.6356
    KmXADXylonate <- 0.7939547312720485
    KmXADkdXyl <- 0.8624275767674491
    KiXADNADH <- 10.462860450236784
    VmKDXD <- 78.857
    KmKDXDkdXyl <- 0.20765619156952056
    KiKDXDxylonate <- 18.30014950694962
    KiKDXDakg <- 14.791937154740701
    KiKDXDpyr <- 17.91684582661257
    KiKDXDlac <- 27.97615617816409
    VmKGSADH <- 58.3
    KmKGSADHGlut <- 0.021750105300746496
    KmKGSADHNAD <- 0.5962856010942978
    KmKGSADHKG <- 0.27905898103297405
    KmKGSADHNADH <- 0.2674
    KmKGSADHKDX <- 0.2136
    KmKGSADHNADH2 <- 0.0241
    KmXLAXylonate <- 0.038100367072300925
    KmKDXDGlut <- 0.28865820908799766
    kLDH <- 10
    fracX <- 0.2511037240696762
    Xyl<-rep(0,NumberOfSteps)
    KG<-rep(0,NumberOfSteps)
    XylLact<-rep(0,NumberOfSteps)
    Xylonate<-rep(0,NumberOfSteps)
    kdXyl<-rep(0,NumberOfSteps)
    Glut<-rep(0,NumberOfSteps)
    NAD<-rep(0,NumberOfSteps)
    pyr<-rep(0,NumberOfSteps)
    lac<-rep(0,NumberOfSteps)
    protb <- InitialGuess[1]
    protc <- InitialGuess[2]
    protd <- InitialGuess[3]
    protx <- InitialGuess[4]
    prota <- InitialGuess[5]
    Ntot<-11.5
    Eb<-protb/1000
    Ec<-protc/1000
    Ed<-protd/1000
    Ex<-0.2511037240696762*protx/1000
    Ea<-prota/1000
    CNE<-1
    Xyl[1]<-4.25
    NAD[1]<-Ntot
    pyr[1]<-10
    Counter<-0
    repeat{
      Counter<-Counter+1

      vXDH<-(VmXDH*Xyl[Counter]*NAD[Counter]/(KmXDHXyl*KmXDHNAD))/((1+Xyl[Counter]/KmXDHXyl+XylLact[Counter]/KmXDHXylLact)*(1+NAD[Counter]/KmXDHNAD+(Ntot-NAD[Counter])/KmXDHNADH))
      vXLANE<-kXLA*XylLact[Counter]
      vXLA<-(VmXLA*XylLact[Counter]/KmXLAXylLact)/(1+XylLact[Counter]/KmXLAXylLact+Xylonate[Counter]/KmXLAXylonate)
      vXAD<-(VmXAD*Xylonate[Counter]/KmXADXylonate)/((1+Xylonate[Counter]/KmXADXylonate+kdXyl[Counter]/KmXADkdXyl)*(1+(Ntot-NAD[Counter])/KiXADNADH))
      vKDXD<-(VmKDXD*(kdXyl[Counter]/KmKDXDkdXyl))/((1+kdXyl[Counter]/KmKDXDkdXyl+Glut[Counter]/KmKDXDGlut)*(1+Xylonate[Counter]/KiKDXDxylonate+KG[Counter]/KiKDXDakg+pyr[Counter]/KiKDXDpyr+lac[Counter]/KiKDXDlac))
      vKGSADH<-(VmKGSADH*Glut[Counter]*NAD[Counter]/(KmKGSADHGlut*KmKGSADHNAD))/((1+Glut[Counter]/KmKGSADHGlut+KG[Counter]/KmKGSADHKG)*(1+NAD[Counter]/KmKGSADHNAD+(Ntot-NAD[Counter])/KmKGSADHNADH)+(kdXyl[Counter]/KmKGSADHKDX)*(1+NAD[Counter]/KmKGSADHNAD+(Ntot-NAD[Counter])/KmKGSADHNADH2))
      vLDH<-kLDH*(Ntot-NAD[Counter])*pyr[Counter]

      dXyl<-(-1*Eb*vXDH)
      dXylLact<-(Eb*vXDH-Ec*vXLA-CNE*vXLANE)
      dXylonate<-(Ec*vXLA+CNE*vXLANE-Ed*vXAD)
      dkdXyl<-(Ed*vXAD-Ex*vKDXD)
      dGlut<-(Ex*vKDXD-Ea*vKGSADH)
      dKG<-(Ea*vKGSADH)
      dNAD<-(-1*Eb*vXDH-Ea*vKGSADH+vLDH)
      dpyr<-(-1*vLDH)
      lac[Counter]<-10-pyr[Counter]
      Xyl[Counter+1]<-Xyl[Counter]+dXyl*dt
      XylLact[Counter+1]<-XylLact[Counter]+dXylLact*dt
      Xylonate[Counter+1]<-Xylonate[Counter]+dXylonate*dt
      kdXyl[Counter+1]<-kdXyl[Counter]+dkdXyl*dt
      Glut[Counter+1]<-Glut[Counter]+dGlut*dt
      KG[Counter+1]<-KG[Counter]+dKG*dt
      NAD[Counter+1]<-NAD[Counter]+dNAD*dt
      pyr[Counter+1]<-pyr[Counter]+dpyr*dt
      if(Counter==NumberOfSteps){break}
    }
    return(-1*KG[NumberOfSteps])
  }

  StartingVector<-c(2,1.5,5.7,9.8,3.6)

  MCGD(5,22.6,50,c(22.6/5,22.6/5,22.6/5,22.6/5,22.6/5),.25,20)
}




