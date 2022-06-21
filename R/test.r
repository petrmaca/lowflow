# library(lowflow)


# dta =read.table("/home/hubert/ownCloud/Shared/petr dp/discharge/daily/Q207290.txt",sep=";",skip=1)

a=0.95
ups =baseflow_UKIH(dta$V4)
ups2=LH_filter(dta$V4,a)
ups3=Chapman_filter(dta$V4,a)

tres = baseflow_RecessionConstant(dta$V4)
tres
bfim=baseflow_BFImax(dta$V4,tres)
ups4=Eckhardt_filter(dta$V4,a=a,BFI_max = bfim)


plot(dta$V4, type="l")
lines(ups,col="red")
lines(ups2,col="blue")
lines(ups3,col="green")
lines(ups3,col="purple")



library(devtools)
install_github("cran/lfstat")
library(lfstat)
wmoBF <- baseflow(dta$V4)
plot(dta$V4, type = "l")
lines(wmoBF, col = 2)

resDF= data.frame(WMO = wmoBF, ukih=ups,LH=ups2,CH =ups3,EK=ups4)
plot(resDF)


resRLM = select.recessionlimbs(dta$V4,minlength = 5)
plot(dta$V4, type = "l")
points(resRLM$t,resRLM$Q,col=2,pch=19)
