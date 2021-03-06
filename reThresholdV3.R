require(tidyr)
require(tibble)
require(readr)
require(ggplot2)
#require(ggrepel)
require(yarrr)
require(pROC)

source("Halo/loadHaloObjectFiles.R")

source("reThresholdFuncs.R")
source("getAllCombinations.R")

args <- commandArgs(trailing=T)
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

pltFile=cc("NewThresholdV3/thresholdROCCurves",
    markerPos,
    paste0(samples,collapse=","),
    ".pdf")

if(!interactive()) pdf(file=pltFile,width=8.5,height=11)

omar=par()$mar

roc.stats.all=list()

if(interactive()) stop("BreakPoint-Alpha")

for(sample in samples) {
    spots=dx %>% filter(Sample==sample) %>% distinct(SPOT) %>% pull(SPOT)
    for(spot in spots) {

        mInten=dd %>%
            filter(Sample==sample & SPOT==spot & Marker==markerPos & ValueType=="Intensity") %>%
            pull(Value)
        amInten=asinh(mInten)

        cat(sample,",",spot,"\n")
        tbl=dd %>%
            filter(Sample==sample & SPOT==spot, ValueType=="Positive") %>%
            count(Marker,Value) %>%
            filter(Value==1) %>%
            arrange(desc(n)) %>%
            filter(Marker %in% c(markerPos,markersToTest))

        roc.stats=list()
        newThetas=c()
        par(mfrow=c(5,2))

        allMarkerNegCombinations=getAllCombinations(markersToTest)
        nCombs=len(allMarkerNegCombinations)
        for(markerNeg in allMarkerNegCombinations) {
            cat("testing",markerNeg,"  ")
            dff=dd %>% filter(Sample==sample & SPOT==spot)
            ss=try({d.roc=getROCMulti(dff,markerPos,markerNeg)})
            if(class(ss)=="roc") {
                stats=getROCStats(sample,spot,markerPos,markerNeg,d.roc)

                newStatI=len(roc.stats)+1
                roc.stats[[newStatI]]=stats

                cat("AUC =",stats$auc,"\n")

                if(!roc.stats[[newStatI]]$thetaOpt %in% newThetas |
                    markerNeg==allMarkerNegCombinations[nCombs]) {
                    newThetas=c(newThetas,roc.stats[[newStatI]]$thetaOpt)
                    par(pty='s')
                    plot.roc.1(d.roc,markerPos,markerNeg)
                    par(pty="m")
                    plot(density(amInten,from=min(amInten),to=max(amInten)),
                        main=paste(markerPos,markerNeg,sample,spot),
                        xlab="asinh(Intensity)")
                    abline(v=asinh(stats$thetaOpt),col="darkgreen",lwd=2,lty=2)
                    abline(v=asinh(thetas[[paste(sample,markerPos,sep=":")]]),col="darkred",lwd=2,lty=2)
                    abline(v=mean(amInten),lty=2,lwd=2,col=8)
                    rug(asinh(d.roc$controls),col="blue",lwd=2)
                    if(markerNeg==allMarkerNegCombinations[nCombs]) {
                        thetaOpt <- roc.stats %>% 
                                        bind_rows %>% 
                                        arrange(desc(auc)) %>% 
                                        slice(1) %>% 
                                        pull(thetaOpt)
                        abline(v=asinh(thetaOpt),col="lightgreen",lwd=2,lty=2)
                    }

                }
            }
        }

        roc.stats.all=c(roc.stats.all,roc.stats)

    }
}

if(!interactive()) dev.off()

roc.tbl=bind_rows(roc.stats.all)
write_csv(roc.tbl,cc("NewThresholdV3/rocStats",paste0(samples,collapse=","),markerPos,".csv"))
write_csv(intensityStats,cc("NewThresholdV3/intensityStats",paste0(samples,collapse=","),markerPos,".csv"))
