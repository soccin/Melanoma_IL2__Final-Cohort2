require(tidyr)
require(tibble)
require(readr)
require(ggplot2)
#require(ggrepel)
require(yarrr)
require(pROC)

source("Halo/loadHaloObjectFiles.R")

source("reThresholdFuncsV4.R")
source("getAllCombinations.R")

markerPos="SOX10"

args <- commandArgs(trailing=T)
doExclusions=TRUE

if(args[1]=="NULL") {
    dblPosMarker=NULL
} else {
    dblPosMarker=args[1]
}
args=args[-1]

i=1
dd=loadHaloObjFile(args[i],exclude=doExclusions)

oDir=file.path("NewThresholdV3",cc(markerPos,dblPosMarker))
dir.create(oDir,showWarnings=F)

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

pltFile=file.path(oDir,cc("thresholdROCCurves",
    markerPos,
    paste0(samples,collapse=","),
    ".pdf"))

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
        theta0=thetas[grep(markerPos,names(thetas))]
        xMax=max(max(mInten),theta0)

        cat(sample,",",spot,"\n")
        tbl=dd %>%
            filter(Sample==sample & SPOT==spot, ValueType=="Positive") %>%
            count(Marker,Value) %>%
            filter(Value==1) %>%
            arrange(desc(n)) %>%
            filter(Marker %in% c(markerPos,markersToTest))

        roc.stats=list()
        newThetas=c()
        if(interactive()) {
            par(mfrow=c(3,2))
        } else {
            par(mfrow=c(5,2))
        }

        allMarkerNegCombinations=getAllCombinations(markersToTest)
        nCombs=len(allMarkerNegCombinations)
        aucMax=0
        for(markerNeg in allMarkerNegCombinations) {
            cat("testing",markerNeg,"  ")
            dff=dd %>% filter(Sample==sample & SPOT==spot)
            numDAPI=dff %>% filter(Marker=="DAPI" & ValueType=="Positive") %>% nrow(.)
            ss=try({d.roc=getROCMulti(dff,markerPos,markerNeg,dblPosMarker)})
            if(class(ss)=="roc") {
                stats=getROCStats(sample,spot,markerPos,markerNeg,d.roc)

                newStatI=len(roc.stats)+1
                roc.stats[[newStatI]]=stats

                cat("AUC =",stats$auc,"\n")

                # if(!roc.stats[[newStatI]]$thetaOpt %in% newThetas |
                #     markerNeg==allMarkerNegCombinations[nCombs]) {
                if(stats$auc>aucMax |
                    markerNeg==allMarkerNegCombinations[nCombs]) {
                    aucMax=stats$auc
                    #newThetas=c(newThetas,roc.stats[[newStatI]]$thetaOpt)
                    par(pty='s')
                    plot.roc.1(d.roc,markerPos,markerNeg)
                    par(pty="m")


                    plot(density(amInten,from=min(amInten),to=max(amInten)),
                        main=paste(markerPos,markerNeg,sample,spot,"\n",
                            paste0("nDAPI =",numDAPI," [SOX10+:",
                            stats$numCases,",negC:",stats$numControls,"]")),
                        xlab="asinh(Intensity)", xlim=c(0,asinh(xMax)))
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
write_csv(roc.tbl,
    file.path(oDir,cc("rocStats",paste0(samples,collapse=","),markerPos,".csv")))
write_csv(intensityStats,
    file.path(oDir,cc("intensityStats",paste0(samples,collapse=","),markerPos,".csv")))

