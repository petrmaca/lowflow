select.recessionlimbs = function(Q,DTM,minlength=9, begOUt = 3, endOUt=5)
{
	reclimbframe = data.frame()
	reclimbnum = 1

	dtm = c()
	tpart = c()
	Qpart = c()
	posinlimbpart = c()

	posinlimbcounter = 1

	N = length(Q)
for(i in 2:(N-1)) {

	if(Q[i+1]<Q[i-1]) {
	
	Qpart = c(Qpart,Q[i])
	tpart = c(tpart,i)

	posinlimbpart = c(posinlimbpart,posinlimbcounter)
	posinlimbcounter = posinlimbcounter+1
} else {
	
	if(length(Qpart)>=minlength) {


	for(j in 1:length(Qpart)) {
		
		reclimbframe = rbind(reclimbframe,
		c(tpart[j],posinlimbpart[j],Qpart[j],reclimbnum))
		}

	reclimbnum = reclimbnum+1
	}

	Qpart = c()
	tpart = c()
	posinlimbpart = c()

	posinlimbcounter = 1
	}
}

names(reclimbframe) = c("t","posinlimb","R","reclimbnum")
reclimbframe = cbind(reclimbframe, as.Date(DTM[reclimbframe$t]))

names(reclimbframe) = c("t","posinlimb","R","reclimbnum","DTM")

nRL = max(reclimbframe$reclimbnum)

reclimbframeCUTED =c()
for( i in 1:nRL){
  helpDF =   reclimbframe[reclimbframe$reclimbnum == i,]
  reclimbframeCUTED = rbind(reclimbframeCUTED,helpDF[(begOUt+1):(nrow(helpDF)-endOUt),])
  rm(helpDF)
}

return(reclimbframeCUTED)

}


select.recessionlimbs_Xie = function(Q,DTM,minlength=9, endOUt=2)
{
  
   # Q = dta1b$R
   # DTM= dta1b$DTM
  reclimbframe = data.frame()
  reclimbnum = 1
  
  dtm = c()
  tpart = c()
  Qpart = c()
  posinlimbpart = c()
  q90 = as.numeric(quantile(dta1b$R_mm_den,probs =c( 0.9)))
  
  posinlimbcounter = 1
  
  N = length(Q)
  for(i in 2:(N-1)) {
    
    if(Q[i+1]<Q[i-1]) {
      
      Qpart = c(Qpart,Q[i])
      tpart = c(tpart,i)
      # dtm = c(dtm,as.character(DTM[i]))
      
      posinlimbpart = c(posinlimbpart,posinlimbcounter)
      posinlimbcounter = posinlimbcounter+1
    } else {
      
      if(length(Qpart)>=minlength) {
        
        
        for(j in 1:length(Qpart)) {
          
          reclimbframe = rbind(reclimbframe,
                               c(dtm[j],tpart[j],posinlimbpart[j],Qpart[j],reclimbnum))
        }
        
        reclimbnum = reclimbnum+1
      }
      
      Qpart = c()
      tpart = c()
      posinlimbpart = c()
      dtm=c()
      
      
      posinlimbcounter = 1
    }
  }
  
  names(reclimbframe) = c("t","posinlimb","R","reclimbnum")
  
  reclimbframe = cbind(reclimbframe, as.Date(DTM[reclimbframe$t]))
  
  names(reclimbframe) = c("t","posinlimb","R","reclimbnum","DTM")
  
  # reclimbframe = cbind(reclimbframe, as.character(dtm))
  
  nRL = max(reclimbframe$reclimbnum)
  
  reclimbframeCUTED =c()
  for( i in 1:nRL){
    helpDF =   reclimbframe[reclimbframe$reclimbnum == i,]
    Qm=max(helpDF$R)
    helpbegOUt <- ifelse(q90<=Qm,5,3)
    reclimbframeCUTED = rbind(reclimbframeCUTED,helpDF[(helpbegOUt+1):(nrow(helpDF)-endOUt),])
    rm(helpDF)
  }
  
  
  return(reclimbframeCUTED)
  
}
