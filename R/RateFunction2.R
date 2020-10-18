#' RateFunction2
#'
#' @param EnzymeConcentration /
#' @param TurnoverRate /
#' @param SubstrateConcentrations /
#' @param MichaelisConstants /
#' @param Substrate1CompetitiveInhibitorConcentrations /
#' @param Substrate1CompetitiveInhibitionConstants /
#' @param Substrate2CompetitiveInhibitorConcentrations /
#' @param Substrate2CompetitiveInhibitionConstants /
#' @param NonCompetitiveInhibitorConcentrations /
#' @param NonCompetitiveInhibitionConstants /

RateFunction2<-function(EnzymeConcentration,TurnoverRate,SubstrateConcentrations,MichaelisConstants,Substrate1CompetitiveInhibitorConcentrations,Substrate1CompetitiveInhibitionConstants,Substrate2CompetitiveInhibitorConcentrations,Substrate2CompetitiveInhibitionConstants,NonCompetitiveInhibitorConcentrations,NonCompetitiveInhibitionConstants){

  Numerator<-TurnoverRate*EnzymeConcentration*SubstrateConcentrations[1]*SubstrateConcentrations[2]

  CIQuotients1<-Substrate1CompetitiveInhibitorConcentrations/Substrate1CompetitiveInhibitionConstants
  CISumTerm1<-sum(CIQuotients1)
  Substrate1Term<-1+SubstrateConcentrations[1]/MichaelisConstants[1]+CISumTerm1

  CIQuotients2<-Substrate2CompetitiveInhibitorConcentrations/Substrate2CompetitiveInhibitionConstants
  CISumTerm2<-sum(CIQuotients2)
  Substrate2Term<-1+SubstrateConcentrations[2]/MichaelisConstants[2]+CISumTerm2

  NonCompetitiveInhibitorQuotients<-NonCompetitiveInhibitorConcentrations/NonCompetitiveInhibitionConstants
  NCISumTerm<-sum(NonCompetitiveInhibitorQuotients)
  NonCompetitiveInhibitionTerm<-1+NCISumTerm

  Denominator<-MichaelisConstants[1]*MichaelisConstants[2]*Substrate1Term*Substrate2Term*NonCompetitiveInhibitionTerm
  Rate<-Numerator/Denominator
  return(Rate)
}
