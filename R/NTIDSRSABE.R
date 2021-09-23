#' NTIDSRSABE
#' An Reference-Scaled Average Bioequivalence Procedure For Therapeutic Index Drug.
#' @param Data is a data.frame,Data must have four column:Subject,Period,Sequence,PK.The contents of the Sequence should be 1 or 2,1 is Sequence of "TRTR",2 is Sequence of "RTRT".
#'
#' @return resultRSABE,the result of RSABE, critbound is upper 90% confidence interval of y.
#' @export NTIDSRSABE
#' @importFrom magrittr %>%
#' @importFrom gmodels estimable
#' @examples
#' NTIDSRSABE(ALLPK)
NTIDSRSABE<-function(Data){
  if ("data.frame" %in% class(Data)) {
    ALLPK = Data
  }
  else {
    stop("Data should be data.frame or file name!")
  }

  ALLPK$Subject<-as.factor(ALLPK$Subject)
  ALLPK$Period<-as.factor(ALLPK$Period)
  ALLPK$Sequence<-as.factor(ALLPK$Sequence)

  SWRtoData<-SWR(ALLPK)
  S2WR<-SWRtoData[1,1]
  DfCVR<-SWRtoData[1,3]

  PKT1_1<-ALLPK %>% filter(Formulation=="T",Sequence==1,Period==1)

  PKT1_2<-ALLPK %>% filter(Formulation=="T",Sequence==2,Period==2)

  PKT1<-rbind(PKT1_1,PKT1_2)
  PKT1<-PKT1[order(PKT1$Subject),]

  PKT2_1<-ALLPK %>% filter(Formulation=="T",Sequence==1,Period==3)

  PKT2_2<-ALLPK %>% filter(Formulation=="T",Sequence==2,Period==4)

  PKT2<-rbind(PKT2_1,PKT2_2)
  PKT2<-PKT2[order(PKT2$Subject),]

  PKR1_1<-ALLPK %>% filter(Formulation=="R",Sequence==1,Period==2)

  PKR1_2<-ALLPK %>% filter(Formulation=="R",Sequence==2,Period==1)

  PKR1<-rbind(PKR1_1,PKR1_2)
  PKR1<-PKR1[order(PKR1$Subject),]

  PKR2_1<-ALLPK %>% filter(Formulation=="R",Sequence==1,Period==4)

  PKR2_2<-ALLPK %>% filter(Formulation=="R",Sequence==2,Period==3)
  PKR2<-rbind(PKR2_1,PKR2_2)
  PKR2<-PKR2[order(PKR2$Subject),]


  PKT1L<-PKT1
  PKT1L$logPKT1<-log(PKT1L$PK)

  PKT2L<-PKT2
  PKT2L$logPKT2<-log(PKT2L$PK)

  PKR1L<-PKR1
  PKR1L$logPKR1<-log(PKR1L$PK)

  PKR2L<-PKR2
  PKR2L$logPKR2<-log(PKR2L$PK)


  PKALLT<-merge(PKT1L,PKT2L,by= "Subject",all=TRUE)

  PKALLR<-merge(PKR1L,PKR2L,by= "Subject",all=TRUE)

  PKALLF<-merge(PKALLT,PKALLR,by = "Subject",all=TRUE)


  PKALLsel<-subset(PKALLF,select = c(Subject,Sequence.x.x,logPKT1,logPKT2,logPKR1,logPKR2))
  names(PKALLsel)[names(PKALLsel) == 'Sequence.x.x'] <- 'Sequence'


  PKALLilat<-PKALLsel

  PKALLilat$ilat<-0.5*(PKALLilat$logPKT1 +PKALLilat$logPKT2 -PKALLilat$logPKR1-PKALLilat$logPKR2)

  PKALLilat$dlat<-(PKALLilat$logPKR1-PKALLilat$logPKR2)

  PKRperiod1<-subset(PKALLilat,logPKR1!="NA")
  PKRperiod2<-subset(PKRperiod1,logPKR2!="NA")
  PKRperiod<-PKRperiod2


  PKTperiod1<-subset(PKALLilat,logPKT1!="NA")
  PKTperiod2<-subset(PKTperiod1,logPKT2!="NA")
  PKTperiod<-PKTperiod2


  PKALLperiod1<-subset(PKALLilat,logPKT1!="NA")
  PKALLperiod2<-subset(PKALLperiod1,logPKT2!="NA")
  PKALLperiod3<-subset(PKALLperiod2,logPKR1!="NA")
  PKALLperiod4<-subset(PKALLperiod3,logPKR2!="NA")
  PKALLperiod<-PKALLperiod4


  PKALLperiod$Subject<-as.factor(PKALLperiod$Subject)
  PKALLperiod$Sequence<-as.factor(PKALLperiod$Sequence)

  modelRilat<- lm(ilat~Sequence, data=PKALLperiod)

  MODEOUT<-gmodels::estimable(modelRilat,cm=c("(Intercept)"=1,Sequence2=0.5),conf.int=0.9)
  Estimate=MODEOUT$Estimate
  SE<-MODEOUT$Std.
  Lower.CI<-MODEOUT$Lower.CI
  Upper.CI<-MODEOUT$Upper.CI
  CIL=exp(MODEOUT$Lower.CI)
  CIU=exp(MODEOUT$Upper.CI)
  Pointest=exp(Estimate)
  x=Estimate**2-SE**2
  boundx<-max(abs(Lower.CI),abs(Upper.CI))**2
  theta=((log(1/0.9))/(0.1))**2
  y=-(theta*S2WR)
  boundy=y*DfCVR/qchisq(0.95,DfCVR)
  critbound=(x+y)+sqrt(((boundx - x)**2) +((boundy - y)**2) )
  names(Estimate) <- c("Estimate")
  names(SE) <- c("SE")
  names(Pointest) <- c("Pointest")
  names(CIL) <- c("CIL")
  names(CIU) <- c("CIU")
  names(x) <- c("x")
  names(boundx) <- c("boundx")
  names(theta) <- c("theta")
  names(y) <- c("y")
  names(boundy) <- c("boundy")
  names(critbound) <- c("critbound")
  resultRSABE<-c(Estimate,SE,Pointest,CIL,CIU,x,boundx,theta,y,boundy,critbound)
  resultRSABE<-data.frame(resultRSABE)
  resultRSABE<-t(resultRSABE)
  return(resultRSABE)
}
