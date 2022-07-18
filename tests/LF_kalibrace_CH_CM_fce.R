rm(list = ls())

require(lowflow)
require(data.table)
require(fst)
require(hydroGOF)
require(DEoptim)

CH = function(Q, a){
  n = length(Q)
  b = vector("numeric", length = n)
  b[1] = Q[1]
  for(i in 2 : n){
    b[i] = b[i - 1] * (3 * a - 1) / (3 - a) + ((1 - a) / (3 - a)) * (Q[i] + Q[i - 1])
    if(b[i] > Q[i]){
      b[i] = Q[i]
      }
    }
  return(b)
  }

CM = function(Q, a){
  n = length(Q)
  b = vector("numeric", length = n)
  b[1] = Q[1]
  for(i in 2 : n){
    b[i] = (a / (2 - a)) * b[i - 1] + ((1 - a) / (2 - a)) * Q[i]
    if(b[i] > Q[i]){
      b[i] = Q[i]
      }
    }
  return(b)
  }

CH_bez_podm = function(Q, a){
  n = length(Q)
  b = vector("numeric", length = n)
  b[1] = Q[1]
  for(i in 2 : n){
    b[i] = b[i - 1] * (3 * a - 1) / (3 - a) + ((1 - a) / (3 - a)) * (Q[i] + Q[i - 1])
    }
  return(b)
  }

CM_bez_podm = function(Q, a){
  n = length(Q)
  b = vector("numeric", length = n)
  b[1] = Q[1]
  for(i in 2 : n){
    b[i] = (a / (2 - a)) * b[i - 1] + ((1 - a) / (2 - a)) * Q[i]
    }
  return(b)
  }

OptimFilter = function(Param_optim){
  alfa = Param_optim[1]
  BFImax = Param_optim[2]
  if(SEL_FIL == 'CH'){
    BF = CH(Q = SEL_POV[, R_mm_den], a = alfa)
    }
  if(SEL_FIL == 'CM'){
    BF = CM(Q = SEL_POV[, R_mm_den], a = alfa)
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

FIL = c('CH', 'CM')
RES_FIL = list()

for(a in 1 : length(FIL)){
  SEL_FIL = FIL[a]

  citac = 1
  RES_PART = data.frame(matrix(nrow = length(IDS), ncol = 8))
  names(RES_PART) = c('FILTR', 'OLA', 'ALFA', 'KGE_podm', 'BIAS_podm', 'KGE_nepodm', 'BIAS_nepodm', 'PROC')
  RES_PART[, 'FILTR'] = SEL_FIL
  
  for(b in 1 : length(IDS)){
    RES_PART[b, 'OLA'] = IDS[b]
    SEL_POV = Q_R[OLA == IDS[b]]
    SEL_POV = SEL_POV[,.(OLA, DTM, R_mm_den)]
    SEL_RECLIMB_ALL = RECLIMBS_obs[OLA == IDS[b]]
    SEL_RECLIMB_ALL = SEL_RECLIMB_ALL[,.(DTM, t, posinlimb, reclimbnum, R)]
    
    if(nrow(SEL_RECLIMB_ALL) > 0){
      
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
        BF = CH(Q = SEL_POV[, R_mm_den], a = alfa)
        }
      if(SEL_FIL == 'CM'){
        BF = CM(Q = SEL_POV[, R_mm_den], a = alfa)
        }
      BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 3))
      names(BF_COMP) = c('OLA', 'DTM', 'selBF')
      BF_COMP[, 'OLA'] = SEL_POV[, OLA]
      BF_COMP[, 'DTM'] = SEL_POV[, DTM]
      BF_COMP[, 'selBF'] = BF
      BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
      RES_PART[b, 'ALFA'] = alfa
      RES_PART[b, 'KGE_podm'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
      RES_PART[b, 'BIAS_podm'] = pbias(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
      
      if(SEL_FIL == 'CH'){
        BF = CH_bez_podm(Q = SEL_POV[, R_mm_den], a = Best_Par[1])
        }
      if(SEL_FIL == 'CM'){
        BF = CM_bez_podm(Q = SEL_POV[, R_mm_den], a = Best_Par[1])
        }
      BF_COMP = data.frame(matrix(nrow = nrow(SEL_POV), ncol = 3))
      names(BF_COMP) = c('OLA', 'DTM', 'selBF')
      BF_COMP[, 'OLA'] = SEL_POV[, OLA]
      BF_COMP[, 'DTM'] = SEL_POV[, DTM]
      BF_COMP[, 'selBF'] = BF
      BF_COMP = merge(BF_COMP, SEL_RECLIMB_ALL, by = 'DTM')
      RES_PART[b, 'KGE_nepodm'] = KGE(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
      RES_PART[b, 'BIAS_nepodm'] = pbias(obs = BF_COMP[, 'R'], sim = BF_COMP[, 'selBF'], na.rm = TRUE)
      BF_COMP[, 'DIF'] = BF_COMP[, 'R'] - BF_COMP[, 'selBF']
      RES_PART[b, 'PROC'] = 100 * length(which(BF_COMP[, 'DIF'] < 0)) / nrow(BF_COMP)
      }
    }
  RES_FIL[[a]] = RES_PART
  }
RES_FIL = rbindlist(RES_FIL)

setwd('C:/prace/LowFlow/Outputs')
saveRDS(RES_FIL, 'kalibrace_CH_CM.rds')

boxplot(RES_FIL[FILTR == 'CH', KGE_podm], RES_FIL[FILTR == 'CM', KGE_podm], outline = FALSE)
boxplot(RES_FIL[FILTR == 'CH', ALFA], RES_FIL[FILTR == 'CM', ALFA], outline = FALSE)
