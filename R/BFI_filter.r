#### functions for BFI #####
## Lyne and Hollick 3passes##
LH_filter_3p= function(Q, a){
  n <- length(Q)
  
  b1 <-vector("numeric", length = n)
  b2 <-vector("numeric", length = n)
  b3 <-vector("numeric", length = n)
  
  f1 <- vector("numeric", length = n)
  f2 <- vector("numeric", length = n)
  f3 <- vector("numeric", length = n)
  
  # the first pass
  f1[1] <- Q[1]
  for(i in 2:n){
    f1[i] <- (a) * f1[i-1] + ((1+a) / (2)) * (Q[i]-Q[i-1])
  }
  f1 <- ifelse(f1 < 0, 0,f1)
  f1 <- ifelse(f1 > Q, Q,f1)
  b1 <- ifelse(f1 > 0, Q-f1,Q )
  
  
  
  # the second pass
  # b2 <- b1
  # f2 <- f1
  f2[n] <- b1[n]
  for(i in (n-1):2){
      f2[i] <- (a) * f2[i+1] + ((1+a) / (2)) * (b1[i]-b1[i+1])
  }
  f2 <- ifelse(f2 < 0, 0,f2)
  f2 <- ifelse(f2 > b1, b1,f2)
  b2 <- ifelse(f2 > 0, b1-f2,b1 )
  

  # the third pass
  f3 <- b2
  for(i in 2:n){
      f3[i] <- (a) * f3[i-1] + ((1+a) / (2)) * (b2[i]-b2[i-1])
  }
  f3 <- ifelse(f3 < 0, 0,f3)
  f3 <- ifelse(f3 > b2, b2,f3)
  b3 <- ifelse(f3 >0, b2-f3,b2 )

  return(b3)
} 
# ekhardt lyne hock 1 pass
LH_filter_1p= function(Q, a){
  n <- length(Q)
  
  b1 <-vector("numeric", length = n)
  
  b1[1] <- Q[1]
  for(i in 2:n){
    b1[i] <- (a) * b1[i-1] + ((1-a) / (2)) * (Q[i]+Q[i-1])
  }
  b1 <- ifelse(b1 >Q, Q,b1 )

  return(b1)
} 



##Chapman Filter which is equal to same filter ith Echartds BFI_MAX = O.5#####
# paper https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/91WR01007
Chapman_filter = function(Q, a){
  n <- length(Q)
  b <-Q
  b[1] <- Q[1]
  f_i_1 <- 0
  for(i in 2:n){
    if(b[i-1] < Q[i]) {
      f<- (3*a-1) / (3-a) * f_i_1 + (2) / (3-a) * (Q[i] - a*Q[i-1])
      b[i] <- (a) * b[i-1] + 1 / 2 *(1-a) * (f + (Q[i-1]-b[i-1]))
      f_i_1 = f
    } else if( b[i-1] > Q[i]){
      b[i] = Q[i]
    }
  }
  return(b)
}

##Chapman MAxwell Filter which is equal to same filter ith Echartds BFI_MAX = O.5#####

Chapman_MAxwell_filter = function(Q, a){
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

obsBaseflow = function(Q,DTM){
  
  
  
  
  
}
