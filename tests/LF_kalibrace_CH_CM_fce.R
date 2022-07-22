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

#vypocet CR a BFImax dle procent
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

#prumerne alfy pro CH a CM
setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
ALFA_CH_CM = readRDS('kalibrace_CH_CM.rds')
PROCENTA = readRDS('Proc_CR_BFImax_BRUT_LANG.rds')
ALFA_PRUM_CH = mean(ALFA_CH_CM[FILTR == 'CH', ALFA])
ALFA_PRUM_CM = mean(ALFA_CH_CM[FILTR == 'CM', ALFA])

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

PROC_BRUT_CH = PRUMERY[which(PRUMERY[, 'DIF_CH_BR'] == min(PRUMERY[, 'DIF_CH_BR'])), 'PROC']
PROC_LANG_CH = PRUMERY[which(PRUMERY[, 'DIF_CH_LN'] == min(PRUMERY[, 'DIF_CH_LN'])), 'PROC']
PROC_BRUT_CM = PRUMERY[which(PRUMERY[, 'DIF_CM_BR'] == min(PRUMERY[, 'DIF_CM_BR'])), 'PROC']
PROC_LANG_CM = PRUMERY[which(PRUMERY[, 'DIF_CM_LN'] == min(PRUMERY[, 'DIF_CM_LN'])), 'PROC']

VEC_K_BRUT_CH = PROCENTA[PROC == PROC_BRUT_CH, RC_Brut]
VEC_K_LANG_CH = PROCENTA[PROC == PROC_LANG_CH, RC_Brut]
VEC_K_BRUT_CM = PROCENTA[PROC == PROC_BRUT_CM, RC_Brut]
VEC_K_LANG_CM = PROCENTA[PROC == PROC_LANG_CM, RC_Brut]

#vypocet lowflow a KGE
RES_KGE = data.frame(matrix(nrow = length(IDS), ncol = 7))
names(RES_KGE) = c('OLA', 'KGE_0925_CH', 'KGE_BRUT_CH', 'KGE_LANG_CH', 'KGE_0925_CM', 'KGE_BRUT_CM', 'KGE_LANG_CM')

for(a in 1 : length(IDS)){
  print(a)
  RES_KGE[a, 'OLA'] = IDS[a]
  SEL_POV = Q_R[OLA == IDS[a]]
  SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den)]
  SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[a]]
  SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
  
  BF_CH = Chapman_filter(Q = SEL_POV[, R_mm_den], a = 0.925)
  BF_BRUT_CH = Chapman_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_BRUT_CH[a])
  BF_LANG_CH = Chapman_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_LANG_CH[a])
  BF_CM = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = 0.925)
  BF_BRUT_CM = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_BRUT_CM[a])
  BF_LANG_CM = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = VEC_K_LANG_CM[a])
  
  BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 8))
  names(BF_COMP) = c('OLA', 'DTM', 'CH', 'BRUT_CH', 'LANG_CH', 'CM', 'BRUT_CM', 'LANG_CM')
  BF_COMP[, 'OLA'] = SEL_POV[, OLA]
  BF_COMP[, 'DTM'] = SEL_POV[, DTM]
  BF_COMP[, 'CH'] = BF_CH
  BF_COMP[, 'BRUT_CH'] = BF_BRUT_CH
  BF_COMP[, 'LANG_CH'] = BF_LANG_CH
  BF_COMP[, 'CM'] = BF_CM
  BF_COMP[, 'BRUT_CM'] = BF_BRUT_CM
  BF_COMP[, 'LANG_CM'] = BF_LANG_CM
  
  BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
  RES_KGE[a, 'KGE_0925_CH'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'CH'], na.rm = TRUE)
  RES_KGE[a, 'KGE_BRUT_CH'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BRUT_CH'], na.rm = TRUE)
  RES_KGE[a, 'KGE_LANG_CH'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'LANG_CH'], na.rm = TRUE)
  RES_KGE[a, 'KGE_0925_CM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'CM'], na.rm = TRUE)
  RES_KGE[a, 'KGE_BRUT_CM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BRUT_CM'], na.rm = TRUE)
  RES_KGE[a, 'KGE_LANG_CM'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'LANG_CM'], na.rm = TRUE)
  }

DTA_graf = data.frame(matrix(nrow = nrow(RES_KGE), ncol = 8))
names(DTA_graf) = c('KGE_0925_CH', 'KGE_kal_CH', 'KGE_BRUT_CH', 'KGE_LANG_CH', 'KGE_0925_CM', 'KGE_kal_CM', 'KGE_BRUT_CM', 'KGE_LANG_CM')
DTA_graf[, 1] = RES_KGE[, 2]
DTA_graf[, 2] = ALFA_CH_CM[FILTR == 'CH', KGE]
DTA_graf[, 3] = RES_KGE[, 3]
DTA_graf[, 4] = RES_KGE[, 4]
DTA_graf[, 5] = RES_KGE[, 5]
DTA_graf[, 6] = ALFA_CH_CM[FILTR == 'CM', KGE]
DTA_graf[, 7] = RES_KGE[, 6]
DTA_graf[, 8] = RES_KGE[, 7]

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
