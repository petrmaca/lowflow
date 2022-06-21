# library(lowflow)
library(data.table)


dtaC = readRDS("tests/data/QRdta151.rds")
dta1b =dtaC[OLA=="UPS003000aqh",]
a=0.95
ups =baseflow_UKIH(dta1b[,R_mm_den])
ups2=LH_filter(dta1b[,R_mm_den],a)
ups3=Chapman_filter(dta1b[,R_mm_den],a)

tres = baseflow_RecessionConstant(dta1b[,R_mm_den])
tres
bfim=baseflow_BFImax(dta1b[,R_mm_den],tres)
ups4=Eckhardt_filter(dta1b[,R_mm_den],a=a,BFI_max = bfim)


plot(dta1b[,R_mm_den], type="l")
lines(ups,col="red")
lines(ups2,col="blue")
lines(ups3,col="green")
lines(ups3,col="purple")



# library(devtools)
# install_github("cran/lfstat")
# library(lfstat)
vmoT=TRUE
if(vmoT){
  wmoBF <- baseflow(dta1b[,R_mm_den])
  plot(dta1b[,R_mm_den], type = "l")
  lines(wmoBF, col = 2)
  
  resDF= data.frame(WMO = wmoBF, ukih=ups,LH=ups2,CH =ups3,EK=ups4)
  plot(resDF)
  
  
  
  
} else {
  
  resDF= data.frame(ukih=ups,LH=ups2,CH =ups3,EK=ups4)
  plot(resDF)
  resRLM = select.recessionlimbs(dta1b[,R_mm_den],minlength = 5)
  plot(dta1b[,R_mm_den], type = "l")
  points(resRLM$t,resRLM$Q,col=2,pch=19)
  
}
