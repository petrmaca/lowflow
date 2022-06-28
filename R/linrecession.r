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
	dtm = c(dtm,as.character(DTM[i]))

	posinlimbpart = c(posinlimbpart,posinlimbcounter)
	posinlimbcounter = posinlimbcounter+1
} else {
	
	if(length(Qpart)>minlength) {


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

names(reclimbframe) = c("DTM","t","posinlimb","Q","reclimbnum")

nRL = max(reclimbframe$reclimbnum)

reclimbframeCUTED =c()
for( i in 1:nRL){
  helpDF =   reclimbframe[reclimbframe$reclimbnum == i,]
  reclimbframeCUTED = rbind(reclimbframeCUTED,helpDF[(begOUt+1):(nrow(helpDF)-endOUt),])
  rm(helpDF)
}

return(reclimbframeCUTED)

}
