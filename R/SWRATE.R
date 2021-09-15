#' SWRATE
#'
#' Evaluate the ratio of SWT/SWR and 90% confidence interval.
#'
#' @param Data is a data.frame,Data must have four column:Subject,Period,Sequence,PK.The contents of the Sequence should be 1 or 2,1 is Sequence of "TRTR",2 is Sequence of "RTRT".
#'
#' @return SWRATE is a data.frame,ResultSWR have show the result of S2WR(The square of the standard deviation within the individual),SWR(standard deviation within the individual),DfCVR(The degree of freedom of standard deviation within an individual),CVwR(Intra-individual variation).
#' @export SWRATE
#' @importFrom magrittr %>%
#'
#'
#' @examples
#' SWRATE(ALLPK)
SWRATE<-function(Data){
  if ("data.frame" %in% class(Data)) {
    ALLPK = Data
  }
  else {
    stop("Data should be data.frame or file name!")
  }

  SWR.Data<-SWR(ALLPK)
  S2WR<-SWR.Data[1,1]
  SWR<-SWR.Data[1,2]
  DFCVR<-SWR.Data[1,3]
  CVwR<-SWR.Data[1,4]

  SWT.Data<-SWT(ALLPK)
  S2WT<-SWT.Data[1,1]
  SWT<-SWT.Data[1,2]
  DFCVT<-SWT.Data[1,3]
  CVwT<-SWT.Data[1,4]
  sw.ratio <- SWT/SWR
  sw.ratio.CI <- c(sw.ratio/sqrt(qf(0.1/2, df1 = DFCVT, df2 = DFCVR,lower.tail = FALSE)),
                   sw.ratio/sqrt(qf(1-0.1/2, df1 = DFCVT, df2 = DFCVR,lower.tail = FALSE)))
  names(sw.ratio.CI) <- c("lower", "upper")
  names(sw.ratio)<-c("sw.ratio")
  names(S2WR)<-c("S2WR")
  names(SWR)<-c("SWR")
  names(DFCVR)<-c("DFCVR")
  names(CVwR)<-c("CVwR")
  names(S2WT)<-c("S2WT")
  names(SWT)<-c("SWT")
  names(DFCVT)<-c("DFCVT")
  names(CVwT)<-c("CVwT")
  SWRATE<-c(S2WR,SWR,DFCVR,CVwR,S2WT,SWT,DFCVT,CVwT,sw.ratio,sw.ratio.CI)
  SWRATE<-data.frame(SWRATE)
  SWRATE<-t(SWRATE)
  return(SWRATE)
}
