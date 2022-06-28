library(lubridate)

# dta1b[,prevDTM := ymd(DTM)-days(1)]
# dta1b[,nextDTM := ymd(DTM)+days(1)]
# helpDT <- dta1b[nextDTM %in% DTM,.(DTM,nextDTM,prevDTM,R_curr = R_mm_den)]
# setkey(helpDT,nextDTM)
# setkey(dta1b,DTM)
# helpdtt<- dta1b[helpDT]
# setkey(helpdtt,i.prevDTM)
# helpDTfin<- dta1b[helpdtt]
# helpDTfin[,dRdt := (i.R_mm_den - R_mm_den)/2]
# finDT<-helpDTfin[,.(DTM = i.DTM.1,nextDTM = nextDTM, prevDTM = prevDTM,Rcur = R_curr,Rnext=i.R_mm_den,Rprev=R_mm_den,dRdt = dRdt) ]
# finDT<-na.omit(finDT)
# finDT[,crit1 := ifelse(dRdt>=0,0,1)]

library(lowflow)
library(data.table)
library(fst)

dtaC = as.data.table(read_fst(path="tests/data/dta.fst"))

vec =dtaC[,unique(ID)]
IDdta <- copy(dtaC[,.(ID = unique(ID),area= unique(AREA))])

IDdta <- IDdta[-c(62,75,92,145,150),]
obsDF =list()
for( i in 1:nrow(IDdta)){
  print(i)
# dt=select.recessionlimbs(finDT$Rcur,as.character(finDT$DTM),minlength=9, begOUt = 3, endOUt=5)
  dta1b =dtaC[ID==IDdta$ID[i],]
  dt=as.data.table(select.recessionlimbs(dta1b$R_mm_den,dta1b$DTM,minlength=9, begOUt = 3, endOUt=5))
  dt[,OLA := unique(dta1b$OLA)]
  dt[,DBCN := unique(dta1b$DBCN)]
  dt[,UPOV_ID := unique(dta1b$UPOV_ID)]
  obsDF[[i]]=dt
}
obsDF <- rbindlist(obsDF)

write_fst(obsDF,path="tests/data/obsRL.fst")
