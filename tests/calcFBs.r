library(lowflow)
library(data.table)
library(fst)
library(hydroGOF)

# 62
# 75
# 92
# 145
# 150
dtaC = as.data.table(read_fst(path="tests/data/dta.fst"))
dtaC[]
vec =dtaC[,unique(ID)]
IDdta <- copy(dtaC[,.(ID = unique(ID),area= unique(AREA))])

IDdta <- IDdta[-c(62,75,92,145,150),]

# copy(dtaC[,.(ID = unique(ID),.(dbcn=DBCN,upovId=UPOV_ID,chp14s=CHP_14_S,xwgs=X_WGS,ywgs =Y_WGS,xutm=X_UTM,yutm=Y_UTM,xjtsk=X_JTSK,yjtsk=Y_JTSK,area=AREA)]
# basinProfiles= copy(dtaC[,.(ID = unique(ID),area=unique(AREA),chp14s=unique(CHP_14_S),xwgs=unique(X_WGS),ywgs =unique(Y_WGS),area= unique(AREA))])
# mybasp =basinProfiles[-c(62,75,92,145,150),]
# saveRDS(mybasp,file="tests/data/petrPovodi.rds")

obsDF  = as.data.table(read_fst(path = "tests/data/obsRL_CenteredDIFFS_3OR5begOUt_2endOut_9length_Xie.fst"))
Id = obsDF[,unique(OLA)]



a = 0.925
nse = c()
for(i in 1:length(Id)){
  # i=1
  basin1D <- obsDF[OLA %in% Id[i],]
  Q1D <- dtaC[OLA %in% Id[i],R_mm_den]
  print(paste(i,nrow(basin1D),Id[i]))
  simBF <- Chapman_filter(Q= Q1D,a)
  bfDF = as.data.table(data.frame(bf = simBF, DTM= dtaC[OLA %in% Id[i],DTM]))
  nse[i] <- NSE(sim=bfDF[DTM %in% basin1D[,DTM],bf], obs =basin1D[,R] )
}
nse
median(nse)


nseVA = c()
avar =c()
for(i in 1:length(Id)){
  # i=1
  basin1D <- obsDF[OLA %in% Id[i],]
  Q1D <- dtaC[OLA %in% Id[i],R_mm_den]
  print(paste(i,nrow(basin1D),Id[i]))
  avar[i] = baseflow_RecessionConstant(Q1D,UB_prc=0.5,method="Langbein")
  simBF <- Chapman_filter(Q= Q1D,avar[i])
  bfDF = as.data.table(data.frame(bf = simBF, DTM= dtaC[OLA %in% Id[i],DTM]))
  nseVA[i] <- NSE(sim=bfDF[DTM %in% basin1D[,DTM],bf], obs =basin1D[,R] )
}
nseVA
mean(avar)
median(nseVA)
median(nse)
mean(nseVA)
mean(nse)


