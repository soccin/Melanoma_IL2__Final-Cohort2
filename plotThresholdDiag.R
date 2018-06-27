
require(tidyverse)
require(fs)

dx=dir_ls("NewThresholdV4/SOX10_S100B/CD20,CD8,CD3,PCK26",reg="rocS.*mel.*csv") %>% map(read_csv) %>% bind_rows
dx=dx %>% mutate(fracPos=numPOS/numDAPI)

ii=dx$MarkerNeg=="CD20:CD8:CD3:PCK26"

pdf(file="testReAssignParams.pdf",width=11,height=6.5)
par(mfrow=c(1,2))

plot(data.frame(dx$fracPos,log10(dx$PCT.Delta)),col=ifelse(ii,"red",8),pch=ifelse(ii,19,1))
abline(v=.02,lty=2,col=8)

plot(dx$auc,dx$precOpt)


with(dx[ii,],points(auc,precOpt,col=2,pch=19,cex=.8))
with(dx[dx$numPOS<500,],points(auc,precOpt,col="darkgreen",pch=19,cex=.8))
with(dx[dx$fracPos<.02,],points(auc,precOpt,col="green",pch=19,cex=.8))
abline(v=.8,lty=2,col=8)
abline(h=.65,lty=2,col=8)


dev.off()
