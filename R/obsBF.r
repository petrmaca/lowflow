library(lubridate)

dta1b[,prevDTM := ymd(DTM)-days(1)]
dta1b[,nextDTM := ymd(DTM)+days(1)]
helpDT <- dta1b[nextDTM %in% DTM,.(DTM,nextDTM,prevDTM,R_curr = R_mm_den)]
setkey(helpDT,nextDTM)
setkey(dta1b,DTM)
helpdtt<- dta1b[helpDT]
setkey(helpdtt,i.prevDTM)
helpDTfin<- dta1b[helpdtt]
helpDTfin[,dRdt := (i.R_mm_den - R_mm_den)/2]
finDT<-helpDTfin[,.(DTM = i.DTM.1,nextDTM = nextDTM, prevDTM = prevDTM,Rcur = R_curr,Rnext=i.R_mm_den,Rprev=R_mm_den,dRdt = dRdt) ]
finDT<-na.omit(finDT)
finDT[,crit1 := ifelse(dRdt>=0,0,1)]


dt=select.recessionlimbs(finDT$Rcur,as.character(finDT$DTM),3)

