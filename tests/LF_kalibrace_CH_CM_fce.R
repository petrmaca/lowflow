rm(list = ls())

require(lowflow)
require(data.table)
require(fst)
require(hydroGOF)
require(DEoptim)

OptimFilter = function(Param_optim){
  alfa = Param_optim[1]
  BFImax = Param_optim[2]
  if(SEL_FIL == 'CH'){
    BF = Chapman_filter(Q = SEL_POV[, R_mm_den], a = alfa)
    }
  if(SEL_FIL == 'CM'){
    BF = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = alfa)
    }
  BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 3))
  names(BF_COMP) = c('OLA', 'DTM', 'selBF')
  BF_COMP[, 'OLA'] = SEL_POV[, OLA]
  BF_COMP[, 'DTM'] = SEL_POV[, DTM]
  BF_COMP[, 'selBF'] = BF
  BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
  CRIT = 1 - KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
  }

Q_R = as.data.table(read_fst(path="tests/data/dta.fst"))
RECLIMBS_obs = as.data.table(read_fst(path="tests/data/obsRL_CenteredDIFFS_3OR5begOUt_2endOut_9length_Xie.fst"))

IDS = unique(Q_R[, OLA])
IDS = IDS[-62]

FIL = c('CH', 'CM')
RES_FIL = list()

#kalibrace alfy pro CH a CM
for(a in 1 : length(FIL)){
  SEL_FIL = FIL[a]

  RES_PART = data.frame(matrix(nrow = length(IDS), ncol = 6))
  names(RES_PART) = c('FILTR', 'OLA', 'ALFA', 'KGE', 'NSE', 'BIAS')
  RES_PART[, 'FILTR'] = SEL_FIL
  
  for(b in 1 : length(IDS)){
    RES_PART[b, 'OLA'] = IDS[b]
    SEL_POV = Q_R[OLA == IDS[b]]
    SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den)]
    SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[b]]
    SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]

    lowerFL = c(0, 0)
    upperFL = c(1, 1)
    np = 20
    pocet_iter = 50
      
    optDE = DEoptim(fn = OptimFilter, 
                    lower = lowerFL, 
                    upper = upperFL, 
                    control = DEoptim.control(VTR = 0, strategy = 4, itermax = pocet_iter, NP = np))
    Best_Par = optDE$optim$bestmem
      
    alfa = Best_Par[1]

    if(SEL_FIL == 'CH'){
      BF = Chapman_filter(Q = SEL_POV[, R_mm_den], a = alfa)
      }
    if(SEL_FIL == 'CM'){
      BF = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = alfa)
      }
    BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 3))
    names(BF_COMP) = c('OLA', 'DTM', 'selBF')
    BF_COMP[, 'OLA'] = SEL_POV[, OLA]
    BF_COMP[, 'DTM'] = SEL_POV[, DTM]
    BF_COMP[, 'selBF'] = BF
    BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
    RES_PART[b, 'ALFA'] = alfa
    RES_PART[b, 'KGE'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
    RES_PART[b, 'NSE'] = NSE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
    RES_PART[b, 'BIAS'] = pbias(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
    }
  RES_FIL[[a]] = RES_PART
  }
RES_FIL = rbindlist(RES_FIL)

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
saveRDS(RES_FIL, 'kalibrace_CH_CM.rds')

#vypocet RC a BFImax dle procent, RC pak pouzity jako alfy pro CH a CM
RES_ALL_POV_PROC = list()

for(a in 1 : length(IDS)){
  print(a)
  SEL_POV = Q_R[OLA == IDS[a]]
  SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den, AREA)]
  
  RES_PART = data.frame(matrix(nrow = 100, ncol = 6))
  names(RES_PART) = c('OLA', 'PROC', 'RC_Brut', 'RC_Lang', 'BFI_max_Brut', 'BFI_max_Lang')
  RES_PART[, 'OLA'] = IDS[a]
  
  for(b in 1 : 100){
    RC_Brut = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = (b / 100), method = 'Brutsaert')
    RC_Lang = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = (b / 100), method = 'Langbein')
    BFI_MAX_BR = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Brut)
    BFI_MAX_LN = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Lang)
    RES_PART[b, 'PROC'] = (b / 100)
    RES_PART[b, 'RC_Brut'] = RC_Brut
    RES_PART[b, 'RC_Lang'] = RC_Lang
    RES_PART[b, 'BFI_max_Brut'] = BFI_MAX_BR
    RES_PART[b, 'BFI_max_Lang'] = BFI_MAX_LN
    }
  RES_ALL_POV_PROC[[a]] = RES_PART
  }
RES_ALL_POV_PROC = rbindlist(RES_ALL_POV_PROC)
setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
saveRDS(RES_ALL_POV_PROC, 'Proc_CR_BFImax_BRUT_LANG.rds')

#vypocet prumerne alfy pro CH a CM z kalibrace
setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
ALFA_CH_CM = readRDS('kalibrace_CH_CM.rds') #alfy z kalibrace
PROCENTA = readRDS('Proc_CR_BFImax_BRUT_LANG.rds') #vypoctene RC a BFImax na zaklade procent
ALFA_PRUM_CH = mean(ALFA_CH_CM[FILTR == 'CH', ALFA])
ALFA_PRUM_CM = mean(ALFA_CH_CM[FILTR == 'CM', ALFA])

#vypocet prumerne alfy pro jednotliva procenta
PRUMERY = data.frame(matrix(nrow = 100, ncol = 7))
names(PRUMERY) = c('PROC', 'PRUM_ALFA_BRUT', 'PRUM_ALFA_LANG', 'DIF_CH_BR', 'DIF_CH_LN', 'DIF_CM_BR', 'DIF_CM_LN')
for(a in 1 : 100){
  PRUMERY[a, 'PROC'] = a / 100
  PRUMERY[a, 'PRUM_ALFA_BRUT'] = mean(PROCENTA[PROC == (a / 100), RC_Brut])
  PRUMERY[a, 'PRUM_ALFA_LANG'] = mean(PROCENTA[PROC == (a / 100), RC_Lang])
  }
PRUMERY[, 'DIF_CH_BR'] = abs(PRUMERY[, 'PRUM_ALFA_BRUT'] - ALFA_PRUM_CH)
PRUMERY[, 'DIF_CH_LN'] = abs(PRUMERY[, 'PRUM_ALFA_LANG'] - ALFA_PRUM_CH)
PRUMERY[, 'DIF_CM_BR'] = abs(PRUMERY[, 'PRUM_ALFA_BRUT'] - ALFA_PRUM_CM)
PRUMERY[, 'DIF_CM_LN'] = abs(PRUMERY[, 'PRUM_ALFA_LANG'] - ALFA_PRUM_CM)

#nalezani procenta, ktere dava prumer alf, ktery je nejpodobnejsi prumeru alf z kalibrace
PROC_BRUT_CH = PRUMERY[which(PRUMERY[, 'DIF_CH_BR'] == min(PRUMERY[, 'DIF_CH_BR'])), 'PROC']
PROC_LANG_CH = PRUMERY[which(PRUMERY[, 'DIF_CH_LN'] == min(PRUMERY[, 'DIF_CH_LN'])), 'PROC']
PROC_BRUT_CM = PRUMERY[which(PRUMERY[, 'DIF_CM_BR'] == min(PRUMERY[, 'DIF_CM_BR'])), 'PROC']
PROC_LANG_CM = PRUMERY[which(PRUMERY[, 'DIF_CM_LN'] == min(PRUMERY[, 'DIF_CM_LN'])), 'PROC']

VEC_K_BRUT_CH = PROCENTA[PROC == PROC_BRUT_CH, RC_Brut]
VEC_K_LANG_CH = PROCENTA[PROC == PROC_LANG_CH, RC_Brut]
VEC_K_BRUT_CM = PROCENTA[PROC == PROC_BRUT_CM, RC_Brut]
VEC_K_LANG_CM = PROCENTA[PROC == PROC_LANG_CM, RC_Brut]

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
saveRDS(PRUMERY, 'Dif_alfa_Brut_Lang.rds')

#nalezeni procenta, ktere dava alfy nejpodobnejsi kalibraci
setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
ALFA_CH_CM = readRDS('kalibrace_CH_CM.rds') #alfy z kalibrace
PROCENTA = readRDS('Proc_CR_BFImax_BRUT_LANG.rds') #vypoctene RC a BFImax na zaklade procent
ALFA_CH_KAL = ALFA_CH_CM[FILTR == 'CH', ALFA]
ALFA_CM_KAL = ALFA_CH_CM[FILTR == 'CM', ALFA]

RES_DIF = data.frame(matrix(nrow = 100, ncol = 5))
names(RES_DIF) = c('PROC', 'CH_Brut', 'CH_Lang', 'CM_Brut', 'CM_Lang')
for(a in 1 : 100){
  RES_DIF[a, 'PROC'] = a / 100
  SEL_PROC = PROCENTA[PROC == (a / 100)]
  RES_DIF[a, 'CH_Brut'] = sum(abs(ALFA_CH_KAL - SEL_PROC[, RC_Brut]))
  RES_DIF[a, 'CH_Lang'] = sum(abs(ALFA_CH_KAL - SEL_PROC[, RC_Lang]))
  RES_DIF[a, 'CM_Brut'] = sum(abs(ALFA_CM_KAL - SEL_PROC[, RC_Brut]))
  RES_DIF[a, 'CM_Lang'] = sum(abs(ALFA_CM_KAL - SEL_PROC[, RC_Lang]))
  }

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
saveRDS(RES_DIF, 'Dif_alfa_Brut_Lang_nejpod_kal.rds')

#vypocet lowflow a KGE pro alfa = 0.925 a alfy urcene jako recesni konstanty dle Brut a Lang 
#(alfy jsou vypocteny na zaklade procenta, ktere vede k alfam, jejichz prumer odpovida prumeru kalibracnich alf, nebo alfy jsou vypocteny na zaklade procenta, ktere dava alfy co nejpodobnejsi kalibracnim alfam)

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
DIF_PRUM = readRDS('Dif_alfa_Brut_Lang.rds')
DIF_NPOD = readRDS('Dif_alfa_Brut_Lang_nejpod_kal.rds')
PROCENTA = readRDS('Proc_CR_BFImax_BRUT_LANG.rds')

PROC_BRUT_CH_PRUM = DIF_PRUM[which(DIF_PRUM[, 'DIF_CH_BR'] == min(DIF_PRUM[, 'DIF_CH_BR'])), 'PROC']
PROC_LANG_CH_PRUM = DIF_PRUM[which(DIF_PRUM[, 'DIF_CH_LN'] == min(DIF_PRUM[, 'DIF_CH_LN'])), 'PROC']
PROC_BRUT_CM_PRUM = DIF_PRUM[which(DIF_PRUM[, 'DIF_CM_BR'] == min(DIF_PRUM[, 'DIF_CM_BR'])), 'PROC']
PROC_LANG_CM_PRUM = DIF_PRUM[which(DIF_PRUM[, 'DIF_CM_LN'] == min(DIF_PRUM[, 'DIF_CM_LN'])), 'PROC']

PROC_BRUT_CH_NPOD = DIF_NPOD[which(DIF_NPOD[, 'CH_Brut'] == min(DIF_NPOD[, 'CH_Brut'])), 'PROC']
PROC_LANG_CH_NPOD = DIF_NPOD[which(DIF_NPOD[, 'CH_Lang'] == min(DIF_NPOD[, 'CH_Lang'])), 'PROC']
PROC_BRUT_CM_NPOD = DIF_NPOD[which(DIF_NPOD[, 'CM_Brut'] == min(DIF_NPOD[, 'CM_Brut'])), 'PROC']
PROC_LANG_CM_NPOD = DIF_NPOD[which(DIF_NPOD[, 'CM_Lang'] == min(DIF_NPOD[, 'CM_Lang'])), 'PROC']

VEC_K_BRUT_CH_PRUM = PROCENTA[PROC == PROC_BRUT_CH_PRUM, RC_Brut]
VEC_K_LANG_CH_PRUM = PROCENTA[PROC == PROC_LANG_CH_PRUM, RC_Brut]
VEC_K_BRUT_CM_PRUM = PROCENTA[PROC == PROC_BRUT_CM_PRUM, RC_Brut]
VEC_K_LANG_CM_PRUM = PROCENTA[PROC == PROC_LANG_CM_PRUM, RC_Brut]

VEC_K_BRUT_CH_NPOD = PROCENTA[PROC == PROC_BRUT_CH_NPOD, RC_Brut]
VEC_K_LANG_CH_NPOD = PROCENTA[PROC == PROC_LANG_CH_NPOD, RC_Brut]
VEC_K_BRUT_CM_NPOD = PROCENTA[PROC == PROC_BRUT_CM_NPOD, RC_Brut]
VEC_K_LANG_CM_NPOD = PROCENTA[PROC == PROC_LANG_CM_NPOD, RC_Brut]

RES_KGE = data.frame(matrix(nrow = length(IDS), ncol = 11))
names(RES_KGE) = c('OLA', 'KGE_0925_CH', 'KGE_BRUT_CH_PRUM', 'KGE_LANG_CH_PRUM', 'KGE_BRUT_CH_NPOD', 'KGE_LANG_CH_NPOD', 'KGE_0925_CM', 'KGE_BRUT_CM_PRUM', 'KGE_LANG_CM_PRUM', 'KGE_BRUT_CM_NPOD', 'KGE_LANG_CM_NPOD')

for(a in 1 : length(IDS)){
  print(a)
  RES_KGE[a, 'OLA'] = IDS[a]
  SEL_POV = Q_R[OLA == IDS[a]]
  SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den)]
  SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[a]]
  SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
  
  BF_CH = Chapman_filter(Q = SEL_POV[, R_mm_den], a = 0.925)
  BF_BRUT_CH_PRUM = Chapman_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_BRUT_CH_PRUM[a])
  BF_LANG_CH_PRUM = Chapman_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_LANG_CH_PRUM[a])
  BF_BRUT_CH_NPOD = Chapman_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_BRUT_CH_NPOD[a])
  BF_LANG_CH_NPOD = Chapman_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_LANG_CH_NPOD[a])
  BF_CM = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = 0.925)
  BF_BRUT_CM_PRUM = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_BRUT_CM_PRUM[a])
  BF_LANG_CM_PRUM = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_LANG_CM_PRUM[a])
  BF_BRUT_CM_NPOD = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_BRUT_CM_NPOD[a])
  BF_LANG_CM_NPOD = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_LANG_CM_NPOD[a])
  
  BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 12))
  names(BF_COMP) = c('OLA', 'DTM', 'CH', 'BRUT_CH_PRUM', 'LANG_CH_PRUM', 'BRUT_CH_NPOD', 'LANG_CH_NPOD', 'CM', 'BRUT_CM_PRUM', 'LANG_CM_PRUM', 'BRUT_CM_NPOD', 'LANG_CM_NPOD')
  BF_COMP[, 'OLA'] = SEL_POV[, OLA]
  BF_COMP[, 'DTM'] = SEL_POV[, DTM]
  BF_COMP[, 'CH'] = BF_CH
  BF_COMP[, 'BRUT_CH_PRUM'] = BF_BRUT_CH_PRUM
  BF_COMP[, 'LANG_CH_PRUM'] = BF_LANG_CH_PRUM
  BF_COMP[, 'BRUT_CH_NPOD'] = BF_BRUT_CH_NPOD
  BF_COMP[, 'LANG_CH_NPOD'] = BF_LANG_CH_NPOD
  BF_COMP[, 'CM'] = BF_CM
  BF_COMP[, 'BRUT_CM_PRUM'] = BF_BRUT_CM_PRUM
  BF_COMP[, 'LANG_CM_PRUM'] = BF_LANG_CM_PRUM
  BF_COMP[, 'BRUT_CM_NPOD'] = BF_BRUT_CM_NPOD
  BF_COMP[, 'LANG_CM_NPOD'] = BF_LANG_CM_NPOD
  
  BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
  RES_KGE[a, 'KGE_0925_CH'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'CH'], na.rm = TRUE)
  RES_KGE[a, 'KGE_BRUT_CH_PRUM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BRUT_CH_PRUM'], na.rm = TRUE)
  RES_KGE[a, 'KGE_LANG_CH_PRUM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'LANG_CH_PRUM'], na.rm = TRUE)
  RES_KGE[a, 'KGE_BRUT_CH_NPOD'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BRUT_CH_NPOD'], na.rm = TRUE)
  RES_KGE[a, 'KGE_LANG_CH_NPOD'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'LANG_CH_NPOD'], na.rm = TRUE)
  RES_KGE[a, 'KGE_0925_CM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'CM'], na.rm = TRUE)
  RES_KGE[a, 'KGE_BRUT_CM_PRUM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BRUT_CM_PRUM'], na.rm = TRUE)
  RES_KGE[a, 'KGE_LANG_CM_PRUM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'LANG_CM_PRUM'], na.rm = TRUE)
  RES_KGE[a, 'KGE_BRUT_CM_NPOD'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BRUT_CM_NPOD'], na.rm = TRUE)
  RES_KGE[a, 'KGE_LANG_CM_NPOD'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'LANG_CM_NPOD'], na.rm = TRUE)
  }

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
saveRDS(RES_KGE, 'KGE_CH_CM_0925_Brut_Lang_proc.rds')

#grafy vysledku
setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
KAL_CH_CM = readRDS('kalibrace_CH_CM.rds')
PROC_CH_CM = readRDS('KGE_CH_CM_0925_Brut_Lang_proc.rds')
DIFS = readRDS('Dif_alfa_Brut_Lang.rds')
PROCENTA = readRDS('Proc_CR_BFImax_BRUT_LANG.rds')

#porovnani alf kal x alf Brut / Lang
ALFA_KAL_CH = KAL_CH_CM[FILTR == 'CH', ALFA]
ALFA_KAL_CM = KAL_CH_CM[FILTR == 'CM', ALFA]

PROC_BRUT_CH = DIFS[which(DIFS[, 'DIF_CH_BR'] == min(DIFS[, 'DIF_CH_BR'])), 'PROC']
PROC_LANG_CH = DIFS[which(DIFS[, 'DIF_CH_LN'] == min(DIFS[, 'DIF_CH_LN'])), 'PROC']
PROC_BRUT_CM = DIFS[which(DIFS[, 'DIF_CM_BR'] == min(DIFS[, 'DIF_CM_BR'])), 'PROC']
PROC_LANG_CM = DIFS[which(DIFS[, 'DIF_CM_LN'] == min(DIFS[, 'DIF_CM_LN'])), 'PROC']

VEC_K_BRUT_CH = PROCENTA[PROC == PROC_BRUT_CH, RC_Brut]
VEC_K_LANG_CH = PROCENTA[PROC == PROC_LANG_CH, RC_Brut]
VEC_K_BRUT_CM = PROCENTA[PROC == PROC_BRUT_CM, RC_Brut]
VEC_K_LANG_CM = PROCENTA[PROC == PROC_LANG_CM, RC_Brut]

plot(ALFA_KAL_CH, VEC_K_BRUT_CH, xlim = c(0.8, 1), ylim = c(0.8, 1), xlab = 'alfa kal.', ylab = 'alfa proc.')
abline(0, 1, col = 'red')
plot(ALFA_KAL_CH, VEC_K_LANG_CH, xlim = c(0.8, 1), ylim = c(0.8, 1), xlab = 'alfa kal.', ylab = 'alfa proc.')
abline(0, 1, col = 'red')
plot(ALFA_KAL_CM, VEC_K_BRUT_CM, xlim = c(0.8, 1), ylim = c(0.8, 1), xlab = 'alfa kal.', ylab = 'alfa proc.')
abline(0, 1, col = 'red')
plot(ALFA_KAL_CM, VEC_K_LANG_CM, xlim = c(0.8, 1), ylim = c(0.8, 1), xlab = 'alfa kal.', ylab = 'alfa proc.')
abline(0, 1, col = 'red')

#porovnani KGE pro al;fa = 0.925, alfa kal. a afly dle procent
DTA_graf = data.frame(matrix(nrow = nrow(PROC_CH_CM), ncol = 8))
names(DTA_graf) = c('0925_CH', 'kal_CH', 'BRUT_CH', 'LANG_CH', '0925_CM', 'kal_CM', 'BRUT_CM', 'LANG_CM')
DTA_graf[, 1] = PROC_CH_CM[, 2]
DTA_graf[, 2] = KAL_CH_CM[FILTR == 'CH', KGE]
DTA_graf[, 3] = PROC_CH_CM[, 3]
DTA_graf[, 4] = PROC_CH_CM[, 4]
DTA_graf[, 5] = PROC_CH_CM[, 5]
DTA_graf[, 6] = KAL_CH_CM[FILTR == 'CM', KGE]
DTA_graf[, 7] = PROC_CH_CM[, 6]
DTA_graf[, 8] = PROC_CH_CM[, 7]

boxplot(DTA_graf, las = 2)

###################



#KGE CH a CM
boxplot(DTA_graf, las = 2)

#porovnani alfy pro CH
plot(ALFA_CH_CM[FILTR == 'CH', ALFA], VEC_K_BRUT_CH, xlab = 'alfa kal.', 
     ylab = 'alfa Brut.', xlim = c(0.85,1), ylim = c(0.85,1))
abline(0, 1, col = 'red')

plot(ALFA_CH_CM[FILTR == 'CH', ALFA], VEC_K_LANG_CH, xlab = 'alfa kal.', 
     ylab = 'alfa Lang', xlim = c(0.85,1), ylim = c(0.85,1))
abline(0, 1, col = 'red')

plot(VEC_K_BRUT_CH, VEC_K_LANG_CH, xlab = 'alfa Brut.', 
     ylab = 'alfa Lang.', xlim = c(0.85,1), ylim = c(0.85,1))
abline(0, 1, col = 'red')

#porovnani alfy pro CM
plot(ALFA_CH_CM[FILTR == 'CM', ALFA], VEC_K_BRUT_CM, xlab = 'alfa kal.', 
     ylab = 'alfa Brut.', xlim = c(0.85,1), ylim = c(0.85,1))
abline(0, 1, col = 'red')

plot(ALFA_CH_CM[FILTR == 'CM', ALFA], VEC_K_LANG_CM, xlab = 'alfa kal.', 
     ylab = 'alfa Lang.', xlim = c(0.85,1), ylim = c(0.85,1))
abline(0, 1, col = 'red')

plot(VEC_K_BRUT_CM, VEC_K_LANG_CM, xlab = 'alfa Brut.', 
     ylab = 'alfa Lang.', xlim = c(0.85,1), ylim = c(0.85,1))
abline(0, 1, col = 'red')
