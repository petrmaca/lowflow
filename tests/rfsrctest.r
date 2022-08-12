library(data.table)
dta <-as.data.table(read.table("tests/data/povodi.csv",sep=",",header=TRUE))

dta
plot(dta[,3:10])
plot(dta[,11:15])
plot(dta[,16:20])
plot(dta[,3:20])
dta[,OLA := paste()]

acka <-readRDS("tests/data/RES/Proc_CR_BFImax_BRUT_LANG.rds")
unique(acka$OLA)

setkey(acka,OLA)
setkey(dta,OLA)
aa <- acka[dta,]
aa[,area:=Shape_Area / 1000/1000]
aa[,delkariverm:=stream_SUM_Shape_Length / 1000]
saveRDS(aa,"tests/data/RES/acka_charpov.rds")

aa[PROC %in% 0.6 ,plot(RC_Brut,Stream_density)]
str(aa)
library(randomForestSRC)
rf <-rfsrc(BFI_max_Brut~slope_MEAN+CN_mean+lesy_podil+Stream_density+area+ornapuda_podil+delkariverm, data = aa[PROC %in% 0.6 ,],ntree=1000)
rf <-rfsrc(RC_Brut~BFI_max_Brut, data = aa[PROC %in% 0.6 ,])
plot(rf)


rfsrc(RC_Brut~., data = aa[PROC %in% 0.6 ,])
