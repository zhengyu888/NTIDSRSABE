#' NTIDSABE
#'An Average Bioequivalence Procedure For Therapeutic Index Drug
#' @param Data is a data.frame,Data must have four column:Subject,Period,Sequence,PK.The contents of the Sequence should be 1 or 2,1 is Sequence of "TRTR",2 is Sequence of "RTRT".
#'
#' @return resultABE,the result of ABE, PointestABE is GMR,CIABE is 90% confidence interval of PointestABE.
#' @export NTIDSABE
#' @importFrom magrittr %>%
#'
#' @examples
#' NTIDSABE(ALLPK)
NTIDSABE<-function(Data){
  if ("data.frame" %in% class(Data)) {
    ALLPK = Data
  }
  else {
    stop("Data should be data.frame or file name!")
  }

  ALLPK$Subject<-as.factor(ALLPK$Subject)
  ALLPK$Period<-as.factor(ALLPK$Period)
  ALLPK$Sequence<-as.factor(ALLPK$Sequence)

  modABE<- lm(log(PK) ~ Sequence + Subject%in%Sequence + Period + Formulation, data=ALLPK)
  EstimateABE<-coef(modABE)[["FormulationT"]]
  PointestABE  <- exp(coef(modABE)[["FormulationT"]])
  CIABE<- as.numeric(exp(confint(modABE, "FormulationT", level=0.9)))
  DFABE<- aov(modABE)[[8]]
  names(EstimateABE) <- c("EstimateABE")
  names(PointestABE) <- c("PointestABE")
  names(CIABE) <- c("lowerABE", "upperABE")
  names(DFABE) <- c("DFABE")
  resultABE<-c(EstimateABE,PointestABE,CIABE,DFABE)
  resultABE<-data.frame(resultABE)
  resultABE<-t(resultABE)
  return(resultABE)
}
