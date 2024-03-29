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

cc=aa[PROC %in% 0.6 ,]
cca=na.omit(cc)
library(randomForestSRC)
rf <-rfsrc(log(RC_Brut)~slope_MEAN+CN_mean+lesy_podil+Stream_density+area+ornapuda_podil+delkariverm, data = cca)
# rf <-rfsrc(RC_Brut~slope_MEAN+CN_mean+lesy_podil+Stream_density+area+ornapuda_podil+delkariverm, data = cca,mtry=4,ntree=10000)
plot(rf)
o <-vimp(rf)
plot(o)
rf


dd <- readRDS("tests/data/RES/kalibrace_CH_CM.rds")
setkey(dd,OLA)
setkey(dta,OLA)
merge(dd,dta)
ee <- dd[dta,]

ee<-na.omit(ee)
ee[,area:=Shape_Area / 1000/1000]
ee[FILTR=="CM",delkariverm:=stream_SUM_Shape_Length / 1000]
ee[,delkariverm:=stream_SUM_Shape_Length / 1000]

rf <-rfsrc(ALFA~slope_STD+elev_STD_1+slope_MEAN+CN_mean+lesy_podil+Stream_density+area+ornapuda_podil+delkariverm, data = ee[FILTR=="CH",],importance = TRUE, splitrule = 'random',ntree=1000,mtry=1)
# rf <-rfsrc(ALFA~slope_STD+elev_STD_1+slope_MEAN+CN_mean+lesy_podil+Stream_density+area+ornapuda_podil+delkariverm, data = ee,importance = TRUE, splitrule = 'random',ntree=200)
so <- subsample(rf)
vimpCI <- extract.subsample(so)$var.jk.sel.Z
vimpCI
plot.subsample(so)
plot(so)

plot(rf)
o <-vimp(rf)
plot(o)
rf

summary(lm(ALFA~, data = ee[FILTR=="CH",]))
o.pred <- predict(object = rf, ee[FILTR=="CH",])
plot(get.mv.predicted(o.pred),ee[FILTR=="CH",ALFA])

cor(c(get.mv.predicted(o.pred)),ee[FILTR=="CH",ALFA])

