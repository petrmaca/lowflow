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

plot(RES_ALL[OLA == IDS[1], RC_Brut] ~ RES_ALL[OLA == IDS[1], PROC], type = 'l', xlab = 'perc.', ylab = 'rec. const.',
     xlim = c(0, 1), ylim = c(0, 1))
for(a in 2 : length(IDS)){
  lines(RES_ALL[OLA == IDS[a], RC_Brut] ~ RES_ALL[OLA == IDS[a], PROC])
  }

plot(RES_ALL[OLA == IDS[1], RC_Lang] ~ RES_ALL[OLA == IDS[1], PROC], type = 'l', xlab = 'perc.', ylab = 'rec. const.',
     xlim = c(0, 1), ylim = c(0, 1))
for(a in 2 : length(IDS)){
  lines(RES_ALL[OLA == IDS[a], RC_Lang] ~ RES_ALL[OLA == IDS[a], PROC])
  }

plot(RES_ALL[OLA == IDS[1], BFI_max_Brut] ~ RES_ALL[OLA == IDS[1], PROC], type = 'l', xlab = 'perc.', ylab = 'BFI max.',
     xlim = c(0, 1), ylim = c(0, 1))
for(a in 2 : length(IDS)){
  lines(RES_ALL[OLA == IDS[a], BFI_max_Brut] ~ RES_ALL[OLA == IDS[a], PROC])
  }
plot(RES_ALL[OLA == IDS[1], BFI_max_Lang] ~ RES_ALL[OLA == IDS[1], PROC], type = 'l', xlab = 'perc.', ylab = 'BFI max.',
     xlim = c(0, 1), ylim = c(0, 1))
for(a in 2 : length(IDS)){
  lines(RES_ALL[OLA == IDS[a], BFI_max_Lang] ~ RES_ALL[OLA == IDS[a], PROC])
  }

PRC = 0.56
plot(RES_ALL[PROC == PRC, RC_Brut],RES_ALL[PROC == PRC, RC_Lang], xlim = c(0.9, 1), ylim = c(0.9, 1), xlab = 'Brutsaert', ylab = 'Langbein')
abline(0, 1, col = 'red')

RES_ALL[PROC == 0.57,mean(RC_Lang)]
LANGens <- RES_ALL[PROC == 0.57,.(OLA,RC_Lang)]

RES_ALL_KGEv = data.frame(matrix(nrow = length(IDS), ncol = 10))
names(RES_ALL_KGEv) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                        'CM_0925', 'EK_080_0925')

RES_ALL_NSEv = data.frame(matrix(nrow = length(IDS), ncol = 10))
names(RES_ALL_NSEv) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                        'CM_0925', 'EK_080_0925')

RES_ALL_PBIASv = data.frame(matrix(nrow = length(IDS), ncol = 10))
names(RES_ALL_PBIASv) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                          'CM_0925', 'EK_080_0925')
for(a in 1 : length(IDS)){
  # a=1
  print(a)
  SEL_POV = Q_R[OLA == IDS[a]]
  SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den, AREA)]
  SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[a]]
  SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
  
  BF = data.frame(matrix(nrow = nrow(SEL_POV), ncol =10))
  names(BF) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                'CM_0925', 'EK_080_0925')
  
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
  BF[, 'LH1_0925'] = LH_filter_1p(Q = SEL_POV[, R_mm_den], a = LANGens[OLA == IDS[a],RC_Lang])
  BF[, 'LH3_0925'] = LH_filter_3p(Q = SEL_POV[, R_mm_den], a = LANGens[OLA == IDS[a],RC_Lang])
  
  # RC_Brut = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = 0.95, method = 'Brutsaert')
  # RC_Lang = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = 0.95, method = 'Langbein')
  # BFI_MAX_BR = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Brut)
  # BFI_MAX_LN = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Lang)
  
  BF[, 'CH_0925'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = LANGens[OLA == IDS[a],RC_Lang])
  # BF[, 'CH_RC_BR_95'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = LANGens[OLA == IDS[a],RC_Lang])
  # BF[, 'CH_RC_LN_95'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang)
  
  BF[, 'CM_0925'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = LANGens[OLA == IDS[a],RC_Lang])
  # BF[, 'CM_RC_BR_95'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut)
  # BF[, 'CM_RC_LN_95'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang)
  
  
  BFI_MAX_LN = baseflow_BFImax(SEL_POV[, R_mm_den], LANGens[OLA == IDS[a],RC_Lang])
  BF[, 'EK_080_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = LANGens[OLA == IDS[a],RC_Lang], BFI_max = BFI_MAX_LN)
  # BF[, 'EK_080_RC_BR_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut, BFI_max = BFI8)
  # BF[, 'EK_080_RC_LN_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang, BFI_max = BFI8)
  
  # BF[, 'EK_BFIBR_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925, BFI_max = BFI_MAX_BR)
  # BF[, 'EK_BFILN_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925, BFI_max = BFI_MAX_LN)
  
  # BF[, 'EK_BFIBR_RC_BR_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut, BFI_max = BFI_MAX_BR)
  # BF[, 'EK_BFILN_RC_LN_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang, BFI_max = BFI_MAX_LN)
  
  BF = merge(BF, SEL_RECLIMB_ALL, by = 'DTM')
  
  RES_ALL_KGEv[a, 'OLA'] = IDS[a]
  RES_ALL_KGEv[a, 'UKIH'] = KGE(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'HYS_f'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'HYS_s'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'HYS_l'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'LH1_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'LH3_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'CH_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CH_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CH_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'CM_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CM_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CM_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'EK_080_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_080_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_080_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFIBR_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFILN_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFIBR_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFILN_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  RES_ALL_NSEv[a, 'OLA'] = IDS[a]
  RES_ALL_NSEv[a, 'UKIH'] = NSE(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'HYS_f'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'HYS_s'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'HYS_l'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'LH1_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'LH3_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'CH_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CH_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CH_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'CM_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CM_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CM_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'EK_080_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_080_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_080_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFIBR_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFILN_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFIBR_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFILN_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  RES_ALL_PBIASv[a, 'OLA'] = IDS[a]
  RES_ALL_PBIASv[a, 'UKIH'] = pbias(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'HYS_f'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'HYS_s'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'HYS_l'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'LH1_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'LH3_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'CH_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CH_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CH_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'CM_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CM_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CM_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'EK_080_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_080_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_080_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFIBR_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFILN_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFIBR_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFILN_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  # RES_ALFA[a, 'OLA'] = IDS[a]
  # RES_ALFA[a, 'ALFA_KONST'] = sel_konst_0925
  # RES_ALFA[a, 'ALFA_BR_95'] = RC_Brut
  # RES_ALFA[a, 'ALFA_LN_95'] = RC_Lang
  # RES_ALFA[a, 'BFImax_BR'] = BFI_MAX_BR
  # RES_ALFA[a, 'BFImax_LN'] = BFI_MAX_LN
  
}

dtaPBIAS =read.table(file="BIAS_obsRECLIMBS_simRECLIMBS.txt",header=TRUE)
RES_ALL[PROC == 0.56,mean(RC_Brut)]

boxplot(RES_ALL_PBIASv[,'EK_080_0925'],dtaPBIAS$EK_080_0925)
boxplot(RES_ALL_PBIASv[,'CM_0925'],dtaPBIAS$CM_0925)
boxplot(RES_ALL_PBIASv[,'CH_0925'],dtaPBIAS$CH_0925)
boxplot(RES_ALL_PBIASv[,'LH1_0925'],dtaPBIAS$LH1_0925)
boxplot(RES_ALL_PBIASv[,'LH3_0925'],dtaPBIAS$LH3_0925)

dtaKGE =read.table(file="KGE_obsRECLIMBS_simRECLIMBS.txt",header=TRUE)
boxplot(RES_ALL_KGEv[,'EK_080_0925'],dtaKGE$EK_080_0925)
boxplot(RES_ALL_KGEv[,'CM_0925'],dtaKGE$CM_0925)
boxplot(RES_ALL_KGEv[,'CH_0925'],dtaKGE$CH_0925)
boxplot(RES_ALL_KGEv[,'LH1_0925'],dtaKGE$LH1_0925)
boxplot(RES_ALL_KGEv[,'LH3_0925'],dtaKGE$LH3_0925)


RES_ALL[PROC == 0.53,mean(RC_Brut)]
Brutens <- RES_ALL[PROC == 0.53,.(OLA,RC_Brut)]

RES_ALL_KGEv = data.frame(matrix(nrow = length(IDS), ncol = 10))
names(RES_ALL_KGEv) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                        'CM_0925', 'EK_080_0925')

RES_ALL_NSEv = data.frame(matrix(nrow = length(IDS), ncol = 10))
names(RES_ALL_NSEv) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                        'CM_0925', 'EK_080_0925')

RES_ALL_PBIASv = data.frame(matrix(nrow = length(IDS), ncol = 10))
names(RES_ALL_PBIASv) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                          'CM_0925', 'EK_080_0925')
for(a in 1 : length(IDS)){
  # a=1
  print(a)
  SEL_POV = Q_R[OLA == IDS[a]]
  SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den, AREA)]
  SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[a]]
  SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
  
  BF = data.frame(matrix(nrow = nrow(SEL_POV), ncol =10))
  names(BF) = c('OLA', 'UKIH', 'HYS_f', 'HYS_s', 'HYS_l', 'LH1_0925', 'LH3_0925', 'CH_0925', 
                'CM_0925', 'EK_080_0925')
  
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
  BF[, 'LH1_0925'] = LH_filter_1p(Q = SEL_POV[, R_mm_den], a = Brutens[OLA == IDS[a],RC_Brut])
  BF[, 'LH3_0925'] = LH_filter_3p(Q = SEL_POV[, R_mm_den], a = Brutens[OLA == IDS[a],RC_Brut])
  
  # RC_Brut = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = 0.95, method = 'Brutsaert')
  # RC_Lang = baseflow_RecessionConstant(Q = SEL_POV[, R_mm_den], UB_prc = 0.95, method = 'Langbein')
  # BFI_MAX_BR = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Brut)
  # BFI_MAX_LN = baseflow_BFImax(SEL_POV[, R_mm_den], RC_Lang)
  
  BF[, 'CH_0925'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = Brutens[OLA == IDS[a],RC_Brut])
  # BF[, 'CH_RC_BR_95'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = LANGens[OLA == IDS[a],RC_Lang])
  # BF[, 'CH_RC_LN_95'] = Chapman_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang)
  
  BF[, 'CM_0925'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = Brutens[OLA == IDS[a],RC_Brut])
  # BF[, 'CM_RC_BR_95'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut)
  # BF[, 'CM_RC_LN_95'] = Chapman_MAxwell_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang)
  
  
  BFI_MAX_LN = baseflow_BFImax(SEL_POV[, R_mm_den], Brutens[OLA == IDS[a],RC_Brut])
  BF[, 'EK_080_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = Brutens[OLA == IDS[a],RC_Brut], BFI_max = BFI_MAX_LN)
  # BF[, 'EK_080_RC_BR_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut, BFI_max = BFI8)
  # BF[, 'EK_080_RC_LN_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang, BFI_max = BFI8)
  
  # BF[, 'EK_BFIBR_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925, BFI_max = BFI_MAX_BR)
  # BF[, 'EK_BFILN_0925'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = sel_konst_0925, BFI_max = BFI_MAX_LN)
  
  # BF[, 'EK_BFIBR_RC_BR_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Brut, BFI_max = BFI_MAX_BR)
  # BF[, 'EK_BFILN_RC_LN_95'] = Eckhardt_filter(Q = SEL_POV[, R_mm_den], a = RC_Lang, BFI_max = BFI_MAX_LN)
  
  BF = merge(BF, SEL_RECLIMB_ALL, by = 'DTM')
  
  RES_ALL_KGEv[a, 'OLA'] = IDS[a]
  RES_ALL_KGEv[a, 'UKIH'] = KGE(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'HYS_f'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'HYS_s'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'HYS_l'] = KGE(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'LH1_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'LH3_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'CH_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CH_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CH_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'CM_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CM_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'CM_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_KGEv[a, 'EK_080_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_080_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_080_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFIBR_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFILN_0925'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFIBR_RC_BR_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_KGE[a, 'EK_BFILN_RC_LN_95'] = KGE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  RES_ALL_NSEv[a, 'OLA'] = IDS[a]
  RES_ALL_NSEv[a, 'UKIH'] = NSE(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'HYS_f'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'HYS_s'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'HYS_l'] = NSE(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'LH1_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'LH3_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'CH_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CH_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CH_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'CM_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CM_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'CM_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_NSEv[a, 'EK_080_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_080_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_080_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFIBR_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFILN_0925'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFIBR_RC_BR_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_NSE[a, 'EK_BFILN_RC_LN_95'] = NSE(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  RES_ALL_PBIASv[a, 'OLA'] = IDS[a]
  RES_ALL_PBIASv[a, 'UKIH'] = pbias(obs = BF[, 'R'], sim = BF[, 'UKIH'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'HYS_f'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_f'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'HYS_s'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_s'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'HYS_l'] = pbias(obs = BF[, 'R'], sim = BF[, 'HYS_l'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'LH1_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'LH1_0925'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'LH3_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'LH3_0925'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'CH_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CH_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CH_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CH_RC_LN_95'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'CM_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CM_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'CM_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'CM_RC_LN_95'], na.rm = TRUE)
  RES_ALL_PBIASv[a, 'EK_080_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_080_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_080_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_080_RC_LN_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFIBR_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFILN_0925'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_0925'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFIBR_RC_BR_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFIBR_RC_BR_95'], na.rm = TRUE)
  # RES_ALL_PBIAS[a, 'EK_BFILN_RC_LN_95'] = pbias(obs = BF[, 'R'], sim = BF[, 'EK_BFILN_RC_LN_95'], na.rm = TRUE)
  
  # RES_ALFA[a, 'OLA'] = IDS[a]
  # RES_ALFA[a, 'ALFA_KONST'] = sel_konst_0925
  # RES_ALFA[a, 'ALFA_BR_95'] = RC_Brut
  # RES_ALFA[a, 'ALFA_LN_95'] = RC_Lang
  # RES_ALFA[a, 'BFImax_BR'] = BFI_MAX_BR
  # RES_ALFA[a, 'BFImax_LN'] = BFI_MAX_LN
  
}

dtaPBIAS =read.table(file="BIAS_obsRECLIMBS_simRECLIMBS.txt",header=TRUE)
# RES_ALL[PROC == 0.56,mean(RC_Brut)]

boxplot(RES_ALL_PBIASv[,'EK_080_0925'],dtaPBIAS$EK_080_0925)
boxplot(RES_ALL_PBIASv[,'CM_0925'],dtaPBIAS$CM_0925)
boxplot(RES_ALL_PBIASv[,'CH_0925'],dtaPBIAS$CH_0925)
boxplot(RES_ALL_PBIASv[,'LH1_0925'],dtaPBIAS$LH1_0925)
boxplot(RES_ALL_PBIASv[,'LH3_0925'],dtaPBIAS$LH3_0925)

dtaKGE =read.table(file="KGE_obsRECLIMBS_simRECLIMBS.txt",header=TRUE)
boxplot(RES_ALL_KGEv[,'EK_080_0925'],dtaKGE$EK_080_0925)
boxplot(RES_ALL_KGEv[,'CM_0925'],dtaKGE$CM_0925)
boxplot(RES_ALL_KGEv[,'CH_0925'],dtaKGE$CH_0925)
boxplot(RES_ALL_KGEv[,'LH1_0925'],dtaKGE$LH1_0925)
boxplot(RES_ALL_KGEv[,'LH3_0925'],dtaKGE$LH3_0925)


library(data.table)
Q_RDT = as.data.table(Q_R)
areas =Q_RDT[OLA %in% IDS,unique(AREA)]
plot(areas,Brutens[,RC_Brut])
points(areas,LANGens[,RC_Lang],col='red')

plot(LANGens[,RC_Lang],Brutens[,RC_Brut])
