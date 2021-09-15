#' SWT.
#'
#' SWT is Evaluate the standard deviation within the individual of the test Drug.
#'
#' @param Data is a data.frame,Data must have four column:Subject,Period,Sequence,PK.The contents of the Sequence should be 1 or 2,1 is Sequence of "TRTR",2 is Sequence of "RTRT".
#'
#' @examples SWT(ALLPK)
#'
#' @return ResultSWT is a data.frame, ResultSWT have show the result of S2WT(The square of the standard deviation within the individual),SWT(standard deviation within the individual),DfCVT(The degree of freedom of standard deviation within an individual),CVwT(Intra-individual variation).
#'
#' @export SWT
#' @importFrom magrittr %>%
#' @importFrom dplyr ungroup
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom PowerTOST mse2CV

SWT<-function(Data){
  if ("data.frame" %in% class(Data)) {
  ALLPK = Data
}
else {
  stop("Data should be data.frame or file name!")
}
  ALLPK$Subject<-as.factor(ALLPK$Subject)
  ALLPK$Period<-as.factor(ALLPK$Period)
  ALLPK$Sequence<-as.factor(ALLPK$Sequence)

  PKT<-ALLPK %>% filter(Formulation=="T")
  PKTdouble <- PKT%>%dplyr::group_by(Subject) %>%dplyr::filter(dplyr::n() == 2) %>%dplyr::ungroup()

  modCVT <- lm(log(PK) ~ Sequence + Subject%in%Sequence + Period,data = PKTdouble)
  aovCVT <- anova(modCVT)
  msewT  <- aovCVT["Residuals", "Mean Sq"]
  S2WT<-msewT
  SWT<-sqrt(S2WT)
  DfCVT  <- aovCVT["Residuals", "Df"]
  CVwT   <- PowerTOST::mse2CV(msewT)
  names(S2WT)<-c("S2WT")
  names(SWT)<-c("SWT")
  names(DfCVT)<-c("DfCVT")
  names(CVwT)<-c("CVwT")
  ResultSWT<-c(S2WT,SWT,DfCVT,CVwT)
  ResultSWT<-data.frame(ResultSWT)
  ResultSWT<-t(ResultSWT)
  return(ResultSWT)
}
