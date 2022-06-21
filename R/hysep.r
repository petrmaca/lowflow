baseflow_HYSEP <- function(Q, area_KM2, method=NULL){
  # R implementation of USGS HYSEP baseflow separation algorithms
  # as described in Pettyjohn & Henning (1979) and implemented 
  # in Sloto & Crouse (1996).
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   area_mi2 = area in square miles
  #   method = HYSEP method; options are "fixed", "sliding", "local"
  #
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  #   https://github.com/samzipper/GlobalBaseflow/blob/master/src/BaseflowSeparationFunctions.R
  ## package dependencies
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  if (sum(is.na(Q))>0){
    stop(paste0(sum(is.na(Q)), " missing data points. You need to gap-fill or remove it."))
  }
  area_mi2 <- 0.386102 * area_KM2
  ## calculate interval width (2N*)
  N <- area_mi2^0.2
  int_width <- 2*floor((2*N)/2)+1  # nearest odd integer to 2N
  if (int_width<3)  int_width <- 3
  if (int_width>11) int_width <- 11
  
  ## calculation depends on method
  if (method=="fixed"){
    
    ## fixed interval of width 2Nstar
    n.ints <- ceiling(length(Q)/int_width)
    ints <- rep(seq(1,n.ints), each=int_width)[1:length(Q)]
    
    # build data frame
    df <- data.frame(int = ints,
                     Q = Q)
    
    # summarize by interval
    df <- df %>% 
      group_by(int) %>% 
      summarize(bf = min(Q)) %>% 
      left_join(df, ., by="int")
    
    return(df$bf)
    
  } else if (method=="sliding"){
    
    ## sliding interval of width 2Nstar
    bf <- rollapply(Q, width=int_width, FUN=min, 
                    align="center", partial=T)
    
    return(bf)
    
  } else if (method=="local"){
    
    ## local minimum
    # local minima are points where Q is equal to the sliding interval minimum
    interval_min <- rollapply(Q, width=int_width, FUN=min,
                              align="center", partial=T)
    i_minima <- which(interval_min==Q)
    
    # interpolate between values using na.approx
    bf <- rep(NaN, length(Q))
    bf[i_minima] <- Q[i_minima]
    if (min(i_minima) != 1) bf[1] <- Q[1]*0.5
    bf <- as.numeric(zoo::na.approx(bf, na.rm=F))
    if (max(i_minima) != length(Q)) bf[(max(i_minima)+1):length(Q)] <- bf[max(i_minima)]
    
    # find any bf>Q and set to Q
    i_tooHigh <- which(bf>Q)
    bf[i_tooHigh] <- Q[i_tooHigh]
    
    return(bf)
    
  } else {
    
    # error
    stop("Wrong or missing method. Choose fixed, sliding, or local")
    
  }
  
}