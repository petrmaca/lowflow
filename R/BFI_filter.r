#### functions for BFI #####
## Lyne and Hollick ##
LH_filter = function(Q, a){
  n <- length(Q)
  b <-Q
  b[1] <- Q[1]
  for(i in 2:n){
    if(b[i-1] < Q[i]) {
      b[i] <- (a) * b[i-1] + ((1-a) / (2)) * (Q[i]+Q[i-1])
    } else {
      b[i]<- Q[i]
    }
  }
  return(b)
} 

##Chapman Filter which is equal to same filter ith Echartds BFI_MAX = O.5#####

Chapman_filter = function(Q, a){
  n <- length(Q)
  b <-Q
  b[1] <- Q[1]
  for(i in 2:n){
    if(b[i-1] < Q[i]) {
      b[i] <- (a) / (2-a) * b[i-1] + (1-a) / (2-a) * (Q[i])
    } else if( b[i-1] > Q[i]){
      b[i] = Q[i]
    }
  }
  return(b)
}

###ECHARDT FILTER ###

Eckhardt_filter = function(Q, a, BFI_max){
  n <- length(Q)
  b <- vector("numeric", length = n)
  b[1] <- Q[1]
  for(i in 2:n){
    if(b[i-1] < Q[i]) {
      b[i] <- ((1-BFI_max)*a*b[i-1] + (1-a)*BFI_max*Q[i])/(1-a*BFI_max)
    } else {
      b[i]<- Q[i]
    }
  }
  return(b)
}


#calculation of BFI max backword filter###### 

# b_filter = function(Q, a){
#   n <- length(Q) 
#   b <- vector("numeric", length = n)
#   b[n] <- Q[n]
#   b[n-1] <- b[n]/a
#   for(j in (n-1):1){
#     if(b[j] < Q[j]){
#       b[j-1] <- (b[j]/a)
#     } else{
#       b[j-1]<- Q[j]
#     }
#   }
#   return(b)
# }

baseflow_UKIH <- function(Q, endrule="NA"){
  # R implementation of UKIH baseflow separation algorithm as
  # described in Piggott et al. (2005)
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   endrule = what to do with endpoints, which will always have NAs?
  #     "NA" (default) = retain NAs
  #     "Q" = use Q of the first/last point
  #     "B" = use bf of the first/last point
  #       
  #
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  # https://github.com/samzipper/GlobalBaseflow/blob/master/src/BaseflowSeparationFunctions.R
  ## package dependencies
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  ## fixed interval of width 5
  int_width <- 5
  n.ints <- ceiling(length(Q)/int_width)
  ints <- rep(seq(1,n.ints), each=int_width)[1:length(Q)]
  
  # build data frame
  df <- data.frame(int = ints,
                   day = seq(1,length(ints)),
                   Q = Q)
  
  # summarize by interval
  df <- 
    df %>% 
    group_by(int) %>% 
    summarize(Qmin = min(Q),
              n_int = sum(is.finite(int))) %>% 
    left_join(df, ., by="int") %>% 
    subset(n_int==int_width)
  
  # extract minimum Qmin for each interval; these are
  # candidates to become turning points
  df.mins <- df[df$Q==df$Qmin, ]
  
  # if there are two minima for an interval (e.g. two 
  # days with same Q), choose the earlier one
  df.mins <- df.mins[!duplicated(df.mins$int), ]
  
  ## determine turning points, defined as:
  #    0.9*Qt < min(Qt-1, Qt+1)
  # do this using a weighted rolling min function
  df.mins$iQmin <- rollapply(df.mins$Qmin, width=3, align="center", 
                             fill=NA, FUN=function(z) which.min(z*c(1,0.9,1)))
  df.mins <- subset(df.mins, is.finite(iQmin))  # get rid of first/last point
  TP.day <- df.mins$day[df.mins$iQmin==2]
  TP.Qmin <- df.mins$Qmin[df.mins$iQmin==2]
  
  if (length(TP.day>1)){
    
    # linearly interpolate to length Q
    bf <- rep(NaN, length(Q))
    bf[TP.day] <- TP.Qmin
    bf <- as.numeric(zoo::na.approx(bf, na.rm=F))
    
    # need to fill in NAs?
    if (endrule=="Q"){
      # start
      bf[1:(TP.day[1]-1)] <- Q[1:(TP.day[1]-1)]
      
      # end
      bf[(TP.day[length(TP.day)]+1):length(Q)] <- 
        Q[(TP.day[length(TP.day)]+1):length(Q)]
      
    } else if (endrule=="B") {
      # start
      bf[1:(TP.day[1]-1)] <- bf[TP.day[1]]
      
      # end
      bf[(TP.day[length(TP.day)]+1):length(Q)] <- 
        bf[TP.day[length(TP.day)]]
      
    } else if (endrule != "NA") {
      
      stop("Invalid endrule")
      
    }
    
  } else {
    bf <- rep(0, length(Q))
  }
  
  # find any bf>Q and set to Q
  i_tooHigh <- which(bf>Q)
  bf[i_tooHigh] <- Q[i_tooHigh]
  return(bf)
}
