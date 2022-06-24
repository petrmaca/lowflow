# library(lowflow)
library(data.table)
library(fst)


# 62
# 75
# 92
# 145
# 150
dtaC = as.data.table(read_fst(path="tests/data/dta.fst"))

vec =dtaC[,unique(ID)]
# vec[62]
# dtaC[ID == as.character(vec[62])]


# dtaC[ID =="158500",]

IDdta <- copy(dtaC[,.(ID = unique(ID),area= unique(AREA))])

IDdta <- IDdta[-c(62,75,92,145,150),]
# dtaC[,.(ID = unique(ID),area= unique(AREA))]
# idd<- "003000"
a=0.925
resDF =list()
for( i in 1:nrow(IDdta)){
print(i)
 i=62
dta1b =dtaC[ID==IDdta$ID[i],]
ups =baseflow_UKIH(dta1b[,R_mm_den])
ups2=LH_filter(dta1b[,R_mm_den],a)
ups3=Chapman_filter(dta1b[,R_mm_den],a)
ups4=Chapman_MAxwell_filter(dta1b[,R_mm_den],a)
tres = baseflow_RecessionConstant(dta1b[,R_mm_den])
# tres
bfim=baseflow_BFImax(dta1b[,R_mm_den],tres)
ups5=Eckhardt_filter(dta1b[,R_mm_den],a=a,BFI_max = bfim)
ups6 =  baseflow_HYSEP(dta1b[,R_mm_den], area_KM2=IDdta[ID==IDdta$ID[i],area], method="fixed")
ups7 =  baseflow_HYSEP(dta1b[,R_mm_den], area_KM2=IDdta[ID==IDdta$ID[i],area], method="sliding")
ups8 =  baseflow_HYSEP(dta1b[,R_mm_den], area_KM2=IDdta[ID==IDdta$ID[i],area], method="local")
 plot(dta1b[,R_mm_den], type="l")
 lines(ups,col="red")
 lines(ups2,col="blue")
# lines(ups3,col="green")
# lines(ups4,col="purple")
# lines(ups5,col="steelblue")
# lines(ups6,col="pink")
# lines(ups7,col="yellow")
# lines(ups8,col="grey")
# resHYSEP =data.frame(fixed = ups6,sliding = ups7, local=ups8)
# plot(resHYSEP)
# # library(devtools)
# # install_github("cran/lfstat")
#  library(lfstat)
# vmoT=TRUE
# if(vmoT){
#   wmoBF <- baseflow(dta1b[,R_mm_den])
#   plot(dta1b[,R_mm_den], type = "l")
#   lines(wmoBF, col = 2)
#   
#   resDF= data.frame(WMO = wmoBF, ukih=ups,LH=ups2,CH =ups3,CM=ups4,E = ups5,fixed =ups6,sliding = ups7, local=ups8)
#   plot(resDF)
#   
#   
# } else {
#   
#   resDF= data.frame(ukih=ups,LH=ups2,CH =ups3,EK=ups4)
#   plot(resDF)
#   resRLM = select.recessionlimbs(dta1b[,R_mm_den],minlength = 5)
#   plot(dta1b[,R_mm_den], type = "l")
#   points(resRLM$t,resRLM$Q,col=2,pch=19)
#   
# }
# resDF= data.frame(WMO = wmoBF, ukih=ups,LH=ups2,CH =ups3,CM=ups4,E = ups5,fixed =ups6,sliding = ups7, local=ups8)
# wmoBF from lfstat similar to UKIH
resDF[[i]]<- data.table(data.frame(ukih=ups,LH=ups2,CH =ups3,CM=ups4,E = ups5,fixed =ups6,sliding = ups7, local=ups8))
}
resDF <- rbindlist(resDF)
cc =cor(resDF,use="pairwise.complete.obs")

heatmap(cc,Rowv = NA, Colv = NA)

# plot(resDF)
# 
library(corrplot)
library(ggplot2)
library(tidyr)
corrplot(cc, method = 'square', order = 'FPC', type = 'lower', diag = FALSE)
corrplot(cc, method = 'square', diag = FALSE, order = 'hclust',
         addrect = 3, rect.col = 'blue', rect.lwd = 3, tl.pos = 'd')


# formatted_cors <- function(df){
#   cors(df) %>%
#     map(~rownames_to_column(.x, var="measure1")) %>%
#     map(~pivot_longer(.x, -measure1, "measure2")) %>% 
#     bind_rows(.id = "id") %>%
#     pivot_wider(names_from = id, values_from = value) %>%
#     mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA)) 
# }
# 
# formatted_cors(mtcars) %>%
#   ggplot(aes(measure1, measure2, col=r)) + ## to get the rect filled
#   geom_tile(col="black", fill="white") +
#   geom_point(aes(size = abs(r)), shape=15) +
#   labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation", title="Correlations in Mtcars") +
#   theme_classic()+
#   scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1))  +
#   scale_x_discrete(expand=c(0,0)) +
#   scale_y_discrete(expand=c(0,0)) +
#   theme(text=element_text(family="Roboto")) +
#   scale_size(range=c(1,11), guide=NULL) 
