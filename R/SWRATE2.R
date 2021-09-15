#' SWRATE2
#' use SWR,DFCVR,SWT,DFCVT Evaluate the ratio of SWT/SWR and 90% confidence interval.
#'
#' @param SWR The standard deviation within the individual of the reference Drug.
#' @param DFCVR The degree of freedom of standard deviation within an individual of the reference Drug.
#' @param SWT The standard deviation within the individual of the test Drug.
#' @param DFCVT The degree of freedom of standard deviation within an individual of the test Drug.
#'
#' @return SWRATE,sw.ratio is the ratio of SWT/SWR,sw.ratio.CI is the 90% confidence interval of sw.ratio.
#' @export SWRATE2
#'
#' @examples
#' SWRATE2(SWR=0.9,DFCVR=0.75,SWT=0.95,DFCVT=0.67)
SWRATE2<-function(SWR,DFCVR,SWT,DFCVT){
  SWT=SWT
  DfCVT=DFCVT
  SWR=SWR
  DfCVR=DFCVR
  sw.ratio <- SWT/SWR
  sw.ratio.CI <- c(sw.ratio/sqrt(qf(0.1/2, df1 = DFCVT, df2 = DFCVR,lower.tail = FALSE)),
                   sw.ratio/sqrt(qf(1-0.1/2, df1 = DFCVT, df2 = DFCVR,lower.tail = FALSE)))
  names(sw.ratio.CI) <- c("lower", "upper")
  names(sw.ratio)<-c("sw.ratio")
  SWRATE<-c(sw.ratio,sw.ratio.CI)
  SWRATE<-data.frame(SWRATE)
  SWRATE<-t(SWRATE)
  return(SWRATE)
}
