require(tidyr)
require(tibble)
require(readr)
require(ggplot2)
#require(ggrepel)
require(yarrr)
require(pROC)

source("Halo/loadHaloObjectFiles.R")

source("reThresholdFuncs.R")

args=commandArgs(trailing=T)
doExclusions=TRUE
i=1
dd=loadHaloObjFile(args[i],exclude=doExclusions)

markerPos="SOX10"

dx=dd %>%
    select(Sample,SPOT,UUID,Marker,ValueType,Value) %>%
    spread(ValueType,Value)

thresholdTbl=dx %>%
    group_by(Sample,Marker) %>%
    summarize(Theta=thresholdIntensity(Intensity,Positive)) %>%
    unite(SampleMarker,Sample,Marker,sep=":",remove=F)

thetas=thresholdTbl$Theta
names(thetas)=thresholdTbl$SampleMarker

dx$Threshold=thetas[paste(dx$Sample,dx$Marker,sep=":")]

intensityStats=dx %>%
    group_by(Sample,Marker) %>%
    summarize(Median=median(Intensity),
        Mean=mean(Intensity),
        Max=max(Intensity),
        MaxN=max(Intensity*(1-Positive)),
        MinP=min(ifelse(Positive,Intensity,Inf)))
intensityStats=full_join(thresholdTbl,intensityStats)

studyName="MelanomaIL2_Final_C2"

samples=dx %>% distinct(Sample) %>% pull(Sample)

markersToTest=scan(cc("reThresholdMarkers_",markerPos,".txt"),"")

stop("XXXX")

pltFile=cc("NewThresholdV2/thresholdROCCurves",
    markerPos,
    paste0(samples,collapse=","),
    ".pdf")

if(!interactive()) pdf(file=pltFile,width=11,height=8.5)

roc.stats=list()

for(sample in samples) {
    spots=dx %>% filter(Sample==sample) %>% distinct(SPOT) %>% pull(SPOT)
    for(spot in spots) {
        cat(sample,",",spot,"\n")
        tbl=dd %>%
            filter(Sample==sample & SPOT==spot, ValueType=="Positive") %>%
            count(Marker,Value) %>%
            filter(Value==1) %>%
            arrange(desc(n)) %>%
            filter(Marker %in% c(markerPos,markersToTest))

        #par(mfrow=c(3,4))
        layout(rbind(c(1,2,3,4),c(5,6,7,7)))
        par(pty='s')

        for(markerNeg in markersToTest) {
            dff=dd %>% filter(Sample==sample & SPOT==spot)
            ss=try({d.roc=getROC(dff,markerPos,markerNeg)})
            if(class(ss)=="roc") {
                theta1=mean(coords(d.roc, "b",ret="t"))
                stats=tibble(
                    Sample=sample,
                    Spot=spot,
                    MarkerPos=markerPos,MarkerNeg=markerNeg,
                    auc=as.numeric(auc(d.roc)),
                    thetaOpt=theta1,
                    specOpt=mean(coords(d.roc,theta1,"t","sp")),
                    sensOpt=mean(coords(d.roc,theta1,"t","se")),
                    precOpt=mean(coords(d.roc,theta1,"t","prec")),
                    accOpt=mean(coords(d.roc,theta1,"t","acc")),
                    theta0=asinh(thetas[paste(sample,markerPos,sep=":")]),
                    thetaOrig=thetas[paste(sample,markerPos,sep=":")],
                    spec0=mean(coords(d.roc,theta0,"t","sp")),
                    sens0=mean(coords(d.roc,theta0,"t","se")),
                    prec0=mean(coords(d.roc,theta0,"t","prec")),
                    acc0=mean(coords(d.roc,theta0,"t","acc")))

                stats$F1=2*stats$precOpt*stats$sensOpt/(stats$precOpt+stats$sensOpt)
                stats$F0=2*stats$prec0*stats$sens0/(stats$prec0+stats$sens0)

                roc.stats[[len(roc.stats)+1]]=stats

                plot.roc.1(d.roc,markerPos,markerNeg)
                if(markerNeg=="CD3") {
                    text(0.1,0,paste0(sample,", ",spot),xpd=T,cex=1.4)
                    statsCD3=stats
                }
            } else {
                plot(1,1,type='n',axes=F,xlab="",ylab="")
                text(1,1,paste("ERROR",markerPos,markerNeg),cex=1.41)
            }
        }

        mInten=dd %>%
            filter(Sample==sample & SPOT==spot & Marker==markerPos & ValueType=="Intensity") %>%
            pull(Value)
        amInten=asinh(mInten)

        par(pty="m")
        plot(density(amInten,from=min(amInten),to=max(amInten)),main=paste(markerPos,"CD3",sample,spot),
            xlab="asinh(Intensity)")
        abline(v=asinh(statsCD3$thetaOpt),col="darkgreen",lwd=2,lty=2)
        abline(v=asinh(thetas[[paste(sample,markerPos,sep=":")]]),col="darkred",lwd=2,lty=2)
        abline(v=mean(amInten),lty=2,lwd=2,col=8)

    }
}

if(!interactive()) dev.off()

roc.tbl=bind_rows(roc.stats)
write_csv(roc.tbl,cc("NewThreshold/rocStats",paste0(samples,collapse=","),markerPos,".csv"))
write_csv(intensityStats,cc("NewThreshold/intensityStats",paste0(samples,collapse=","),markerPos,".csv"))
