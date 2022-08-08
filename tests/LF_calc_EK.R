rm(list = ls())

require(lowflow)
require(data.table)
require(fst)
require(hydroGOF)

EK_PODM_cit = function(Q, a, BFImax){
  n = length(Q)
  b = matrix(nrow = n, ncol = 2)
  b[, 2] = 0
  b[1, 1] = Q[1]
  for(i in 2 : n){
    b[i, 1] = ((1 - BFImax) * a * b[i - 1, 1] + (1 - a) * BFImax * Q[i]) / (1 - a * BFImax)
    if(b[i, 1] > Q[i]){
      b[i, 2] = 1
      b[i, 1] = Q[i]
      }
    }
  return(b)
  }

EK_bezPodm = function(Q, a, BFImax){
  n = length(Q)
  b = vector("numeric", length = n)
  b[1] = Q[1]
  for(i in 2 : n){
    b[i] = ((1 - BFImax) * a * b[i - 1] + (1 - a) * BFImax * Q[i]) / (1 - a * BFImax)
    }
  return(b)
  }


Q_R = as.data.table(read_fst(path="tests/data/dta.fst"))
RECLIMBS_obs = as.data.table(read_fst(path="tests/data/obsRL_CenteredDIFFS_3OR5begOUt_2endOut_9length_Xie.fst"))

IDS = unique(Q_R[, OLA])
IDS = IDS[-62]

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
PROCENTA = readRDS('Proc_CR_BFImax_BRUT_LANG.rds')

RES_ALL = list()

for(a in 1 : 100){
  print(a)
  VEC_BRUT_ALFA = PROCENTA[PROC == (a / 100), RC_Brut]
  VEC_LANG_ALFA = PROCENTA[PROC == (a / 100), RC_Lang]
  VEC_BRUT_BFIm = PROCENTA[PROC == (a / 100), BFI_max_Brut]
  VEC_LANG_BFIm = PROCENTA[PROC == (a / 100), BFI_max_Lang]
  
  RES_PART = data.frame(matrix(nrow = length(IDS), ncol = 6))
  names(RES_PART) = c('OLA', 'PROC', 'KGE_Brut', 'KGE_Lang', 'Proc_podm_Brut', 'Proc_podm_Lang')
  RES_PART[, 'PROC'] = (a / 100)
  
  for(b in 1 : length(IDS)){
    RES_PART[b, 'OLA'] = IDS[b]
    SEL_POV = Q_R[OLA == IDS[b]]
    SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den)]
    SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[b]]
    SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
    
    BF_BRUT = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = VEC_BRUT_ALFA[b], BFI_max = VEC_BRUT_BFIm[b])
    BF_LANG = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = VEC_LANG_ALFA[b], BFI_max = VEC_LANG_BFIm[b])
    BF_PODM_BRUT = EK_PODM_cit(Q = SEL_POV[, R_mm_den], a = VEC_BRUT_ALFA[b], BFImax = VEC_BRUT_BFIm[b])
    BF_PODM_LANG = EK_PODM_cit(Q = SEL_POV[, R_mm_den], a = VEC_LANG_ALFA[b], BFImax = VEC_LANG_BFIm[b])
    BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 8))
    names(BF_COMP) = c('OLA', 'DTM', 'BF_BRUT', 'BF_LANG', 'BF_BRUT_podm', 'BF_LANG_podm', 'DIF_Brut', 'DIF_Lang')
    BF_COMP[, 'OLA'] = SEL_POV[, OLA]
    BF_COMP[, 'DTM'] = SEL_POV[, DTM]
    BF_COMP[, 'BF_BRUT'] = BF_BRUT
    BF_COMP[, 'BF_LANG'] = BF_LANG
    BF_COMP[, 'BF_BRUT_podm'] = BF_PODM_BRUT
    BF_COMP[, 'BF_LANG_podm'] = BF_PODM_LANG
    BF_COMP[, 'DIF_Brut'] = SEL_POV[, R_mm_den] - BF_COMP[, 'BF_BRUT_podm']
    BF_COMP[, 'DIF_Lang'] = SEL_POV[, R_mm_den] - BF_COMP[, 'BF_LANG_podm']
    BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
    RES_PART[b, 'KGE_Brut'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BF_BRUT'], na.rm = TRUE)
    RES_PART[b, 'KGE_Lang'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BF_LANG'], na.rm = TRUE)
    RES_PART[b, 'Proc_podm_Brut'] = 100 * length(which(BF_COMP[, 'DIF_Brut'] < 0)) / nrow(BF_COMP)
    RES_PART[b, 'Proc_podm_Lang'] = 100 * length(which(BF_COMP[, 'DIF_Lang'] < 0)) / nrow(BF_COMP)
    }
  RES_ALL[[a]] = RES_PART
  }
RES_ALL = rbindlist(RES_ALL)

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
saveRDS(RES_ALL, 'KGE_EK_Brut_Lang.rds')

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
INP_DTA = readRDS('KGE_EK_Brut_Lang.rds')
PROCENTA = readRDS('Proc_CR_BFImax_BRUT_LANG.rds')
MED_KGE_PRUM_ALFA_BFIm = data.frame(matrix(nrow = 100, ncol = 9))
names(MED_KGE_PRUM_ALFA_BFIm) = c('PROC', 'medKGE_Brut', 'medKGE_Lang', 'prumALFA_Brut', 'prumBFIm_Brut', 'proc_Brut', 'prumALFA_Lang', 'prumBFIm_Lang', 'proc_Lang')

b = 1
PROCT = 0.01
SEL_POV = Q_R[OLA == IDS[b]]
SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den)]
BF_nPOD = EK_bezPodm(Q = SEL_POV[, R_mm_den], a = PROCENTA[OLA == IDS[b] & PROC == PROCT, RC_Brut], BFImax = PROCENTA[OLA == IDS[b] & PROC == PROCT, BFI_max_Brut])
BF_POD = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = PROCENTA[OLA == IDS[b] & PROC == PROCT, RC_Brut], BFI_max = PROCENTA[OLA == IDS[b] & PROC == PROCT, BFI_max_Brut])
print(PROCENTA[OLA == IDS[b] & PROC == PROCT, BFI_max_Brut])
print(PROCENTA[OLA == IDS[b] & PROC == PROCT, RC_Brut])
rng = 1 : 200
plot(SEL_POV[rng, R_mm_den], type = 'l', col = 'red')
lines(BF_nPOD[rng], col = 'blue')
lines(BF_POD[rng], col = 'green')

for(a in 1 : 100){
  MED_KGE_PRUM_ALFA_BFIm[a, 'PROC'] = a / 100
  MED_KGE_PRUM_ALFA_BFIm[a, 'medKGE_Brut'] = median(INP_DTA[PROC == (a / 100), KGE_Brut])
  MED_KGE_PRUM_ALFA_BFIm[a, 'medKGE_Lang'] = median(INP_DTA[PROC == (a / 100), KGE_Lang])
  MED_KGE_PRUM_ALFA_BFIm[a, 'prumALFA_Brut'] = mean(PROCENTA[PROC == (a / 100), RC_Brut])
  MED_KGE_PRUM_ALFA_BFIm[a, 'prumBFIm_Brut'] = mean(PROCENTA[PROC == (a / 100), BFI_max_Brut])
  MED_KGE_PRUM_ALFA_BFIm[a, 'proc_Brut'] = mean(INP_DTA[PROC == (a / 100), Proc_podm_Brut])
  MED_KGE_PRUM_ALFA_BFIm[a, 'prumALFA_Lang'] = mean(PROCENTA[PROC == (a / 100), RC_Lang])
  MED_KGE_PRUM_ALFA_BFIm[a, 'prumBFIm_Lang'] = mean(PROCENTA[PROC == (a / 100), BFI_max_Lang])
  MED_KGE_PRUM_ALFA_BFIm[a, 'proc_Lang'] = mean(INP_DTA[PROC == (a / 100), Proc_podm_Lang])
  }

VEC_BRUT_ALFA = PROCENTA[PROC == 0.53, RC_Brut]
VEC_LANG_ALFA = PROCENTA[PROC == 0.57, RC_Lang]
VEC_BRUT_BFIm = PROCENTA[PROC == 0.53, BFI_max_Brut]
VEC_LANG_BFIm = PROCENTA[PROC == 0.57, BFI_max_Lang]

RES_KGE = data.frame(matrix(nrow = length(IDS), ncol = 4))
names(RES_KGE) = c('OLA', 'KGE_konst', 'KGE_Brut', 'KGE_Lang')

for(a in 1 : length(IDS)){
  
  SEL_POV = Q_R[OLA == IDS[a]]
  SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den)]
  SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[a]]
  SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
  
  BF_KONST = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = 0.925, BFI_max = 0.8)
  BF_BRUT_VAR = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = VEC_BRUT_ALFA[a], BFI_max = VEC_BRUT_BFIm[a])
  BF_LANG_VAR = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = VEC_LANG_ALFA[a], BFI_max = VEC_LANG_BFIm[a])
  
  BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 5))
  names(BF_COMP) = c('OLA', 'DTM', 'BF_KONST', 'BF_BRUT_VAR', 'BF_LANG_VAR')
  BF_COMP[, 'OLA'] = SEL_POV[, OLA]
  BF_COMP[, 'DTM'] = SEL_POV[, DTM]
  BF_COMP[, 'BF_KONST'] = BF_KONST
  BF_COMP[, 'BF_BRUT_VAR'] = BF_BRUT_VAR
  BF_COMP[, 'BF_LANG_VAR'] = BF_LANG_VAR
  BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
  RES_KGE[a, 'OLA'] = IDS[a]
  RES_KGE[a, 'KGE_konst'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BF_KONST'], na.rm = TRUE)
  RES_KGE[a, 'KGE_Brut'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BF_BRUT_VAR'], na.rm = TRUE)
  RES_KGE[a, 'KGE_Lang'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'BF_LANG_VAR'], na.rm = TRUE)
  }



boxplot(RES_KGE[, (2:4)], outline = FALSE)

plot(MED_KGE_PRUM_ALFA_BFIm[, 'medKGE_Brut'] ~ MED_KGE_PRUM_ALFA_BFIm[, 'PROC'], type = 'l', xlab = 'proc.', ylab = 'KGE', col = 'blue')
lines(MED_KGE_PRUM_ALFA_BFIm[, 'medKGE_Lang'] ~ MED_KGE_PRUM_ALFA_BFIm[, 'PROC'], col = 'red')

plot(MED_KGE_PRUM_ALFA_BFIm[, 'prumALFA_Brut'] ~ MED_KGE_PRUM_ALFA_BFIm[, 'PROC'],  type = 'l', xlab = 'proc.', ylab = 'alfa', col = 'red', ylim = c(0.5, 1))
lines(MED_KGE_PRUM_ALFA_BFIm[, 'prumALFA_Lang'] ~ MED_KGE_PRUM_ALFA_BFIm[, 'PROC'], col = 'blue')

plot(MED_KGE_PRUM_ALFA_BFIm[, 'prumBFIm_Brut'] ~ MED_KGE_PRUM_ALFA_BFIm[, 'PROC'],  type = 'l', xlab = 'proc.', ylab = 'BFImax', col = 'red', ylim = c(0, 1))
lines(MED_KGE_PRUM_ALFA_BFIm[, 'prumBFIm_Lang'] ~ MED_KGE_PRUM_ALFA_BFIm[, 'PROC'], col = 'blue')


