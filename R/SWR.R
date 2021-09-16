#' SWR
#'SWR is Evaluate the standard deviation within the individual of the reference Drug.
#' @param Data is a data.frame,Data must have four column:Subject,Period,Sequence,PK.The contents of the Sequence should be 1 or 2,1 is Sequence of "TRTR",2 is Sequence of "RTRT".
#'
#' @return ResultSWR is a data.frame, ResultSWR have show the result of S2WR(The square of the standard deviation within the individual),SWR(standard deviation within the individual),DfCVR(The degree of freedom of standard deviation within an individual),CVwR(Intra-individual variation).
#' @export SWR
#' @importFrom magrittr %>%
#' @importFrom dplyr ungroup
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr filter
#' @importFrom PowerTOST mse2CV
#'
#' @examples
#' SWR(ALLPK)
SWR<-function(Data){
  if ("data.frame" %in% class(Data)) {
    ALLPK = Data
  }
  else {
    stop("Data should be data.frame or file name!")
  }
  ALLPK$Subject<-as.factor(ALLPK$Subject)
  ALLPK$Period<-as.factor(ALLPK$Period)
  ALLPK$Sequence<-as.factor(ALLPK$Sequence)

  PKR<-ALLPK %>% filter(Formulation=="R")
  PKRdouble <- PKR%>%dplyr::group_by(Subject) %>%dplyr::filter(dplyr::n() == 2) %>%dplyr::ungroup()

  modCVR <- lm(log(PK) ~ Sequence + Subject%in%Sequence + Period,data = PKRdouble)
  aovCVR <- anova(modCVR)
  msewR  <- aovCVR["Residuals", "Mean Sq"]
  S2WR<-msewR
  SWR<-sqrt(S2WR)
  DFCVR  <- aovCVR["Residuals", "Df"]
  CVwR   <- PowerTOST::mse2CV(msewR)
  names(S2WR)<-c("S2WR")
  names(SWR)<-c("SWR")
  names(DFCVR)<-c("DFCVR")
  names(CVwR)<-c("CVwR")
  ResultSWR<-c(S2WR,SWR,DFCVR,CVwR)
  ResultSWR<-data.frame(ResultSWR)
  ResultSWR<-t(ResultSWR)
  return(ResultSWR)
}
