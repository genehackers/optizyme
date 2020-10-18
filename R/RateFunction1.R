#' RateFunction1
#'
#' @param EnzymeConcentration /
#' @param TurnoverRate /
#' @param SubstrateConcentration /
#' @param MichaelisConstant /
#' @param CompetitiveInhibitorConcentrations /
#' @param CompetitiveInhibitorConstants /
#' @param NonCompetitiveInhibitorConcentrations /
#' @param NonCompetitiveInhibitionConstants /
#' @return Rate

RateFunction1<-function(EnzymeConcentration,TurnoverRate,SubstrateConcentration,MichaelisConstant,CompetitiveInhibitorConcentrations,CompetitiveInhibitionConstants,NonCompetitiveInhibitorConcentrations,NonCompetitiveInhibitionConstants){

  Numerator<-TurnoverRate*EnzymeConcentration*SubstrateConcentration

  CompetitiveInhibitorQuotients<-CompetitiveInhibitorConcentrations/CompetitiveInhibitionConstants
  CISumTerm<-sum(CompetitiveInhibitorQuotients)
  CompetitiveInhibitionTerm<-1+SubstrateConcentration/MichaelisConstant+CISumTerm

  NonCompetitiveInhibitorQuotients<-NonCompetitiveInhibitorConcentrations/NonCompetitiveInhibitionConstants
  NCISumTerm<-sum(NonCompetitiveInhibitorQuotients)
  NonCompetitiveInhibitionTerm<-1+NCISumTerm

  Denominator<-MichaelisConstant*CompetitiveInhibitionTerm*NonCompetitiveInhibitionTerm
  Rate<-Numerator/Denominator
  return(Rate)

}
