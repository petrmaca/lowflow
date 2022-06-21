baseflow_BFImax <- function(Q, k){
  # Estimate BFImax parameter for Eckhardt baseflow separation filter
  # using a backwards-looking filter, based on Collischonn & Fan (2013).
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   k = recession constant; this can be estimated with the function baseflow_RecessionConstant.
  #
  # Outputs:
  #   BFImax = maximum allowed value of baseflow index; Eckhardt estimates values of:
  #      0.8 for perennial stream with porous aquifer
  #      0.5 for ephemeral stream with porous aquifer
  #      0.25 for perennial stream with hardrock aquifer
  #    based on a few streams in eastern US
  # https://github.com/samzipper/GlobalBaseflow/blob/master/src/BaseflowSeparationFunctions.R
  
  # start from end of timeseries
  bf <- rep(NaN, length(Q))
  bf[length(Q)] <- Q[length(Q)]
  for (i in (length(Q)-1):1){
    if (bf[i+1]==0){
      bf[i] <- Q[i]
    } else {
      bf[i] <- bf[i+1]/k
    }
    
    # make sure bf isn't > Q
    if (bf[i]>Q[i]) bf[i] <- Q[i]
  }
  
  BFImax <- sum(bf)/sum(Q)
  return(BFImax)
  
}

## example data
run.example <- F
if (run.example){
  # packages required for sample data/examples
  require(dataRetrieval)
  require(lubridate)
  require(ggplot2)
  require(magrittr)
  require(reshape2)
  require(EcoHydRology)
  
  # get USGS data for sample site
  dv <- readNWISdv(siteNumber="04148000",                          # site code (can be a vector of multiple sites)
                   parameterCd="00060",                            # parameter code: "00060" is cubic ft/sec
                   startDate="1900-01-01",endDate="2000-12-31",    # start & end dates (YYYY-MM-DD format)
                   statCd = "00003")                               # statistic code: "00003" is daily mean (default)
  colnames(dv) <- c("agency_cd", "site_no", "Date", "discharge.cfs", "QA.code")
  area_mi2 <- 593
  
  # check for missing data
  sum(is.na(dv$discharge.cfs))
  
  # estimate recession constant
  k <- baseflow_RecessionConstant(dv$discharge.cfs, method="Langbein")
  BFImax <- baseflow_BFImax(Q=dv$discharge.cfs, k=k)
  
  ## perform baseflow separations
  dv$HYSEP_fixed <- baseflow_HYSEP(Q = dv$discharge.cfs, area_mi2 = area_mi2, method="fixed")
  dv$HYSEP_slide <- baseflow_HYSEP(Q = dv$discharge.cfs, area_mi2 = area_mi2, method="sliding")
  dv$HYSEP_local <- baseflow_HYSEP(Q = dv$discharge.cfs, area_mi2 = area_mi2, method="local")
  dv$UKIH <- baseflow_UKIH(Q = dv$discharge.cfs, endrule="B")
  dv$BFLOW_1pass <- baseflow_BFLOW(Q = dv$discharge.cfs, beta=0.925, passes=1)
  dv$BFLOW_3pass <- baseflow_BFLOW(Q = dv$discharge.cfs, beta=0.925, passes=3)
  dv$Eckhardt <- baseflow_Eckhardt(Q = dv$discharge.cfs, BFImax=BFImax, k=k)
  
  dv.melt <- 
    dv %>% 
    subset(select=c("Date", "discharge.cfs", "HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH", "BFLOW_1pass", "BFLOW_3pass", "Eckhardt")) %>% 
    melt(id=c("Date", "discharge.cfs"))
  
  ## calculate BFI
  dv.melt %>% 
    group_by(variable) %>% 
    summarize(discharge.sum = sum(discharge.cfs),
              baseflow.sum = sum(value),
              BFI = round(baseflow.sum/discharge.sum, 2))
  
  ## estimates for 04148000 from Neff et al (2005) and Eckhardt (2008)
  #' HYSEP_fixed = 0.71
  #' HYSEP_slide = 0.70
  #' HYSEP_local = 0.57  # slight difference- I think this is due to the treatment of multiple consecutive minimum values that are equal to each other.
  #' UKIH        = 0.53
  #' BFLOW_1pass = 0.70
  #' BFLOW_3pass = 0.47
  #' Eckhardt    = 0.69  # difference because Eckhardt uses BFImax=0.8
  
  ## make plot
  p.date.start <- ymd("1945-01-01")
  p.date.end <- ymd("1946-01-01")
  
  p <- 
    ggplot(subset(dv.melt, Date >= p.date.start & Date <= p.date.end)) +
    geom_ribbon(data=subset(dv, Date >= p.date.start & Date <= p.date.end), 
                aes(x=Date, ymin=0, ymax=discharge.cfs), fill="black") +
    geom_line(aes(x=Date, y=value, color=variable)) +
    scale_y_continuous(name="Discharge [cfs]") +
    scale_x_date(expand=c(0,0)) +
    scale_color_discrete(name="Method") +
    theme_bw() + 
    theme(panel.grid=element_blank(),
          legend.position=c(0.99, 0.99),
          legend.justification=c(1,1))
  p
  
}

if (run.example){
  ## Lochsa River example from Sanchez-Murillo et al paper
  # packages required for sample data/examples
  require(dataRetrieval)
  require(lubridate)
  require(ggplot2)
  require(magrittr)
  require(reshape2)
  
  # get USGS data for sample site
  dv <- readNWISdv(siteNumber="13337000",                          # site code (can be a vector of multiple sites)
                   parameterCd="00060",                            # parameter code: "00060" is cubic ft/sec
                   startDate="1939-01-01",endDate="2011-12-31",    # start & end dates (YYYY-MM-DD format)
                   statCd = "00003")                               # statistic code: "00003" is daily mean (default)
  colnames(dv) <- c("agency_cd", "site_no", "Date", "discharge.cfs", "QA.code")
  area_mi2 <- 1178
  area_km2 <- area_mi2/(0.621371*0.621371)
  
  # calculate Q in mm/d
  cfs_to_mmd <- ((86400/(5280*5280*5280))/area_mi2)*1.60934*1000*1000
  Q <- dv$discharge.cfs*cfs_to_mmd
  
  # calculate dQ/dt
  # calculate difference
  dQ_dt = c(NaN, diff(Q))
  
  # find days of five consecutive negative values
  which_negative <- which(dQ_dt < 0 & Q > 0)
  which_positive <- which(dQ_dt >= 0)
  which_positive_with_buffer <- unique(c(#which_positive-2, which_positive-1,
    which_positive, 
    which_positive+1, which_positive+2, which_positive+3))  # 2 days before and 3 days after a positive or 0 value
  which_positive_with_buffer <- which_positive_with_buffer[which_positive_with_buffer > 0]  # get rid of negative indices
  which_keep <- which_negative[!(which_negative %in% which_positive_with_buffer)]
  
  # make plot in Fig 1
  data.frame(Q = Q[which_keep], 
             dQ_dt = dQ_dt[which_keep]) %>% 
    ggplot(aes(x=Q, y=-dQ_dt)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10()
  
  # estimate recession constant
  k <- baseflow_RecessionConstant(dv$discharge.cfs, method="Langbein")
  BFImax <- baseflow_BFImax(Q=dv$discharge.cfs, k=k)
  
}

if (run.example){
  ## Smith Creek example from Miller et al (2016) WRR
  # packages required for sample data/examples
  require(dataRetrieval)
  require(lubridate)
  require(ggplot2)
  require(magrittr)
  require(reshape2)
  require(dplyr)
  
  # get USGS data for sample site
  dv <- readNWISdv(siteNumber="01632900",                          # site code (can be a vector of multiple sites)
                   parameterCd="00060",                            # parameter code: "00060" is cubic ft/sec
                   startDate="1960-01-01",endDate="2013-09-30",    # start & end dates (YYYY-MM-DD format)
                   statCd = "00003")                               # statistic code: "00003" is daily mean (default)
  colnames(dv) <- c("agency_cd", "site_no", "Date", "discharge.cfs", "QA.code")
  area_mi2 <- 93.6
  
  # calculate Q in mm/d
  cfs_to_mmd <- ((86400/(5280*5280*5280))/area_mi2)*1.60934*1000*1000
  Q <- dv$discharge.cfs*cfs_to_mmd
  
  # calculate dQ/dt
  # calculate difference
  dQ_dt = c(NaN, diff(Q))
  
  # find days of five consecutive negative values
  which_negative <- which(dQ_dt < 0 & Q > 0)
  which_positive <- which(dQ_dt >= 0)
  which_positive_with_buffer <- unique(c(which_positive-2, which_positive-1,
                                         which_positive, 
                                         which_positive+1, which_positive+2))  # 2 days before and 3 days after a positive or 0 value
  which_positive_with_buffer <- which_positive_with_buffer[which_positive_with_buffer > 0]  # get rid of negative indices
  which_keep <- which_negative[!(which_negative %in% which_positive_with_buffer)]
  
  # plot data
  data.frame(Q = Q[which_keep], 
             dQ_dt = dQ_dt[which_keep]) %>% 
    ggplot(aes(x=Q, y=-dQ_dt)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10()
  
  ## results from paper:
  # k = 0.980
  # BFImax = 0.632
  # mean annual BFI = 0.780
  #   I THINK THIS REPORTED BFI VALUE IS WRONG - my results (~0.580) agree exactly with 
  #   the Purdue WHAT results run using Miller's k and BFImax values
  #
  # results from WHAT: 843635.55/1454109.53 = 0.580
  k_brut <- baseflow_RecessionConstant(Q[which_keep], UB_prc=0.90, method="Brutsaert")
  k_lang <- baseflow_RecessionConstant(Q[which_keep], method="Langbein")
  
  k_brut
  k_lang
  
  BFImax <- baseflow_BFImax(Q, k_brut)
  BFImax
  
  eckhardt <- baseflow_Eckhardt(Q, BFImax, k_brut)
  eckhardt <- baseflow_Eckhardt(Q, 0.632, 0.980)
  sum(eckhardt)/sum(Q)
  
  dv.yr <- 
    data.frame(Date = dv$Date,
               year = year(dv$Date),
               Q = Q,
               B = eckhardt) %>% 
    group_by(year) %>% 
    summarize(Qyr = sum(Q),
              Byr = sum(B),
              BFIyr = sum(Byr/Qyr))
  mean(dv.yr$BFIyr)
  
}


baseflow_RecessionConstant <- function(Q, UB_prc=0.95, method="Brutsaert", min_pairs=50){
  # Script to estimate baseflow recession constant.
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   UB_prc = percentile to use for upper bound of regression
  #   method = method to use to calculate recession coefficient
  #     "Langbein" = Langbein (1938) as described in Eckhardt (2008)
  #     "Brutsaert" = Brutsaert (2008) WRR
  #   min_pairs = minimum number of date pairs retained after filtering out 
  #     quickflow events; 50 is from van Dijk (2010) HESS
  #       
  # Output:
  #   k = recession constant
  # https://github.com/samzipper/GlobalBaseflow/blob/master/src/BaseflowSeparationFunctions.R
  ## package dependencies
  require(quantreg)  # used for quantile regression
  
  if (method=="Langbein"){
    # calculate difference
    dQ_dt = c(NaN, diff(Q))
    
    # find days of five consecutive negative values
    which_negative <- which(dQ_dt < 0 & Q > 0)
    which_positive <- which(dQ_dt >= 0)
    which_positive_with_buffer <- unique(c(which_positive-2, which_positive-1,
                                           which_positive, 
                                           which_positive+1, which_positive+2))  # 2 days before and 2 days after a positive or 0 value
    which_positive_with_buffer <- which_positive_with_buffer[which_positive_with_buffer > 0]  # get rid of negative indices
    which_keep <- which_negative[!(which_negative %in% which_positive_with_buffer)]  # get rid of points within buffer around flow increases
    which_keep <- which_keep[(which_keep-1) %in% which_keep]  # trim to dates with both the current and previous day retained
    
    # any data exist to fit?
    if (length(which_keep) >= min_pairs){
      
      # fit regression
      fit.qr <- rq(Q[which_keep] ~ 0+Q[which_keep-1], tau=UB_prc)  # force intercept to go through origin
      
      # extract constant
      k <- as.numeric(coef(fit.qr)[1])
      
    } else {
      k <- NaN
    }
    return(k)
  }
  
  if (method=="Brutsaert"){
    # calculate lagged difference (dQ/dt) based on before/after point
    dQ_dt <- c(NaN, diff(Q, lag=2)/2, NaN)
    dQ_dt_left <- c(NaN, diff(Q))
    
    # screen data for which dQ_dt to calculate recession, based on rules in Brutsaert (2008) WRR Section 3.2
    which_negative <- which(dQ_dt < 0 & dQ_dt_left < 0 & Q > 0)
    which_positive <- which(dQ_dt >= 0)
    which_positive_with_buffer <- unique(c(which_positive-2, which_positive-1, which_positive,
                                           which_positive+1, which_positive+2, which_positive+3))  # 2 days before and 3 days after a positive or 0 value
    which_positive_with_buffer <- which_positive_with_buffer[which_positive_with_buffer > 0]  # get rid of negative indices; possible because of 2 days before
    which_keep <- which_negative[!(which_negative %in% which_positive_with_buffer)]  # get rid of points within buffer around flow increases
    which_keep <- which_keep[(which_keep-1) %in% which_keep]  # trim to dates with both the current and previous day retained
    
    # any data exist to fit?
    if (length(which_keep) >= min_pairs){
      
      # fit regression
      fit.qr <- rq(Q[which_keep] ~ 0+Q[which_keep-1], tau=UB_prc)  # force intercept to go through origin
      
      # extract constant
      k <- as.numeric(coef(fit.qr)[1])
      
    } else {
      k <- NaN
    }
    return(k)
  }
  
  
}