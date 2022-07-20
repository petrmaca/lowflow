rm(list = ls())

require(lowflow)
require(data.table)
require(fst)
require(hydroGOF)

Q_R = as.data.table(read_fst(path="tests/data/dta.fst"))
RECLIMBS_obs = as.data.table(read_fst(path="tests/data/obsRL_CenteredDIFFS_3OR5begOUt_2endOut_9length_Xie.fst"))

IDS = unique(Q_R[, OLA])
IDS = IDS[-62]

RES_ALL = list()

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
  RES_ALL[[a]] = RES_PART
  }
RES_ALL = rbindlist(RES_ALL)

setwd('C:/Users/hermanovsky/Documents/dev/lowflow/tests/data/RES')
saveRDS(RES_ALL, 'Dif_Proc_RC_BFImax_All_Cat.rds')
