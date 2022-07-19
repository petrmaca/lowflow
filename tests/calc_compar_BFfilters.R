rm(list = ls())

require(lowflow)
require(data.table)
require(fst)
require(hydroGOF)

Q_R = as.data.table(read_fst(path="tests/data/dta.fst"))
RECLIMBS_obs = as.data.table(read_fst(path="tests/data/obsRL_CenteredDIFFS_3OR5begOUt_2endOut_9length_Xie.fst"))

IDS = unique(Q_R[, OLA])
IDS = IDS[-62]

RES_ALFA = data.frame(matrix(nrow = length(IDS), ncol = 6))
names(RES_ALFA) = c('OLA', 'ALFA_KONST', 'ALFA_BR_95', 'ALFA_LN_95', 'BFImax_BR', 'BFImax_LN')

RES_ALL_KGE = data.frame(matrix(nrow = length(IDS), ncol = 20))
names(RES_ALL_KGE) = c('OLA','UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                       'CH_RC_BR_95', 'CH_RC_LN_95', 'CM_0925', 'CM_RC_BR_95', 'CM_RC_LN_95',
                       'EK_080_0925', 'EK_080_RC_BR_95', 'EK_080_RC_LN_95', 'EK_BFIBR_0925', 'EK_BFILN_0925',
                       'EK_BFIBR_RC_BR_95', 'EK_BFILN_RC_LN_95')

RES_ALL_NSE = data.frame(matrix(nrow = length(IDS), ncol = 20))
names(RES_ALL_KGE) = c('OLA','UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                       'CH_RC_BR_95', 'CH_RC_LN_95', 'CM_0925', 'CM_RC_BR_95', 'CM_RC_LN_95',
                       'EK_080_0925', 'EK_080_RC_BR_95', 'EK_080_RC_LN_95', 'EK_BFIBR_0925', 'EK_BFILN_0925',
                       'EK_BFIBR_RC_BR_95', 'EK_BFILN_RC_LN_95')

RES_ALL_PBIAS = data.frame(matrix(nrow = length(IDS), ncol = 20))
names(RES_ALL_KGE) = c('OLA','UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                       'CH_RC_BR_95', 'CH_RC_LN_95', 'CM_0925', 'CM_RC_BR_95', 'CM_RC_LN_95',
                       'EK_080_0925', 'EK_080_RC_BR_95', 'EK_080_RC_LN_95', 'EK_BFIBR_0925', 'EK_BFILN_0925',
                       'EK_BFIBR_RC_BR_95', 'EK_BFILN_RC_LN_95')

sel_konst_0925 = 0.925
BFI8 = 0.8

for(a in 1 : length(IDS)){
  print(a)
  SEL_POV = Q_R[OLA == IDS[a]]
  SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den, AREA)]
  SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[a]]
  SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
  
  BF = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 22))
  names(BF) = c('OLA', 'DTM', 'R_mm_den','UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                'CH_RC_BR_95', 'CH_RC_LN_95', 'CM_0925', 'CM_RC_BR_95', 'CM_RC_LN_95', 'EK_080_0925',
                'EK_080_RC_BR_95', 'EK_080_RC_LN_95', 'EK_BFIBR_0925', 'EK_BFILN_0925', 'EK_BFIBR_RC_BR_95',
                'EK_BFILN_RC_LN_95')
  
  BF[, 'OLA'] = IDS[a]
  BF[, 'DTM'] = SEL_POV[, DTM]
  BF[, 'R_mm_den'] = SEL_POV[, R_mm_den]
  #UKIH
  BF[, 'UKIH'] = baseflow_UKIH(Q = SEL_POV[, R_mm_den])
  #graficke metody
  BF[, 'HYS_f'] = baseflow_HYSEP(Q = SEL_POV[, R_mm_den], area_KM2 = unique(SEL_POV[, AREA]), method = 'fixed')
  BF[, 'HYS_s'] = baseflow_HYSEP(Q = SEL_POV[, R_mm_den], area_KM2 = unique(SEL_POV[, AREA]), method = 'sliding')
  BF[, 'HYS_l'] = baseflow_HYSEP(Q = SEL_POV[, R_mm_den], area_KM2 = unique(SEL_POV[, AREA]), method = 'local')
  #digitalni filtry
  BF[, 'LH1_0925'] = LH_filter_1p(Q = SEL_POV[, R_mm_den], a = sel_konst_0925)
  BF[, 'LH3_0925'] = LH_filter_3p(Q = SEL_POV[, R_mm_den], a = sel_konst_0925)
  
  RC_Brut = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = 0.95, method = 'Brutsaert')
  RC_Lang = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = 0.95, method = 'Langbein')
  BFI_MAX_BR = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Brut)
  BFI_MAX_LN = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Lang)
  
  BF[, 'CH_0925'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925)
  BF[, 'CH_RC_BR_95'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut)
  BF[, 'CH_RC_LN_95'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang)
  
  BF[, 'CM_0925'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925)
  BF[, 'CM_RC_BR_95'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut)
  BF[, 'CM_RC_LN_95'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang)
  
  BF[, 'EK_080_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925, BFI_max = BFI8)
  BF[, 'EK_080_RC_BR_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut, BFI_max = BFI8)
  BF[, 'EK_080_RC_LN_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang, BFI_max = BFI8)
  
  BF[, 'EK_BFIBR_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925, BFI_max = BFI_MAX_BR)
  BF[, 'EK_BFILN_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925, BFI_max = BFI_MAX_LN)
  
  BF[, 'EK_BFIBR_RC_BR_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut, BFI_max = BFI_MAX_BR)
  BF[, 'EK_BFILN_RC_LN_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang, BFI_max = BFI_MAX_LN)
  
  BF = merge(BF, SEL_RECLIMB_ALL, by = 'DTM')
  
  RES_ALL_KGE[a, 'OLA'] = IDS[a]
  RES_ALL_KGE[a, 'UKIH'] = KGE(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_KGE[a, 'HYS_f'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_KGE[a, 'HYS_s'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_KGE[a, 'HYS_l'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_KGE[a, 'LH1_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_KGE[a, 'LH3_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_KGE[a, 'CH_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  RES_ALL_KGE[a, 'CH_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  RES_ALL_KGE[a, 'CH_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_KGE[a, 'CM_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  RES_ALL_KGE[a, 'CM_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  RES_ALL_KGE[a, 'CM_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_KGE[a, 'EK_080_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  RES_ALL_KGE[a, 'EK_080_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  RES_ALL_KGE[a, 'EK_080_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  RES_ALL_KGE[a, 'EK_BFIBR_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  RES_ALL_KGE[a, 'EK_BFILN_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  RES_ALL_KGE[a, 'EK_BFIBR_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  RES_ALL_KGE[a, 'EK_BFILN_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  RES_ALL_NSE[a, 'OLA'] = IDS[a]
  RES_ALL_NSE[a, 'UKIH'] = NSE(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_NSE[a, 'HYS_f'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_NSE[a, 'HYS_s'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_NSE[a, 'HYS_l'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_NSE[a, 'LH1_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_NSE[a, 'LH3_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_NSE[a, 'CH_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  RES_ALL_NSE[a, 'CH_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  RES_ALL_NSE[a, 'CH_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_NSE[a, 'CM_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  RES_ALL_NSE[a, 'CM_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  RES_ALL_NSE[a, 'CM_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_NSE[a, 'EK_080_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  RES_ALL_NSE[a, 'EK_080_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  RES_ALL_NSE[a, 'EK_080_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  RES_ALL_NSE[a, 'EK_BFIBR_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  RES_ALL_NSE[a, 'EK_BFILN_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  RES_ALL_NSE[a, 'EK_BFIBR_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  RES_ALL_NSE[a, 'EK_BFILN_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  RES_ALL_PBIAS[a, 'OLA'] = IDS[a]
  RES_ALL_PBIAS[a, 'UKIH'] = pbias(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'HYS_f'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'HYS_s'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'HYS_l'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'LH1_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'LH3_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'CH_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'CH_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'CH_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'CM_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'CM_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'CM_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'EK_080_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'EK_080_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'EK_080_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'EK_BFIBR_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'EK_BFILN_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'EK_BFIBR_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  RES_ALL_PBIAS[a, 'EK_BFILN_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  RES_ALFA[a, 'OLA'] = IDS[a]
  RES_ALFA[a, 'ALFA_KONST'] = sel_konst_0925
  RES_ALFA[a, 'ALFA_BR_95'] = RC_Brut
  RES_ALFA[a, 'ALFA_LN_95'] = RC_Lang
  RES_ALFA[a, 'BFImax_BR'] = BFI_MAX_BR
  RES_ALFA[a, 'BFImax_LN'] = BFI_MAX_LN
  
  }

par(mar = c(9, 6, 2, 2))
boxplot(RES_ALL_KGE[, 2:20], outline = FALSE, las = 2)
boxplot(RES_ALFA[, 2:6], outline = FALSE, las = 2)
