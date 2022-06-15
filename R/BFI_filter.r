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


