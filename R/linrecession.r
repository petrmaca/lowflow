select.recessionlimbs = function(Q,minlength=3)
{
	reclimbframe = data.frame()
	reclimbnum = 1

	tpart = c()
	Qpart = c()
	posinlimbpart = c()

	posinlimbcounter = 1

	N = length(Q)
for(i in 2:N) {

	if(Q[i]<Q[i-1]) {
	
	Qpart = c(Qpart,Q[i-1])
	tpart = c(tpart,i-1)

	posinlimbpart = c(posinlimbpart,posinlimbcounter)
	posinlimbcounter = posinlimbcounter+1
} else {
	
	if(length(Qpart)>minlength) {


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

names(reclimbframe) = c("t","posinlimb","Q","reclimbnum")

return(reclimbframe)

}