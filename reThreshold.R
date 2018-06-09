require(tidyr)
require(tibble)
require(readr)
require(ggplot2)
#require(ggrepel)
require(yarrr)
require(pROC)

source("Halo/loadHaloObjectFiles.R")

getROC<-function(df,markerPos,markerNeg) {

    dt=df %>%
        filter(Marker %in% c(markerPos,markerNeg)) %>%
        unite(MarkerValue,Marker,ValueType) %>%
        select(Sample,SPOT,UUID,MarkerValue,Value) %>%
        spread(MarkerValue,Value)

    posIntensity=cc(markerPos,"Intensity")
    posPositive=cc(markerPos,"Positive")
    negIntensity=cc(markerNeg,"Intensity")
    negPositive=cc(markerNeg,"Positive")

    dt[["posIntensity"]]=dt[[posIntensity]]
    dt[["posPositive"]]=dt[[posPositive]]
    dt[["negIntensity"]]=dt[[negIntensity]]
    dt[["negPositive"]]=dt[[negPositive]]

    dt %<>%
        mutate(category=case_when(
                negPositive==1 ~ "NEG",
                TRUE ~ "POS"))

    d.roc=roc(dt$category,dt$posIntensity)
}

plot.roc.1<-function(d.roc,markerPos,markerNeg) {

    plot(d.roc,main=paste0(markerPos,", ",markerNeg),asp='s',
        print.thres=T,print.thres.cex=1.12,
        print.auc=T,print.auc.cex=1.12)

}

thresholdIntensity<-function(ii,p){sum(c(min(ii[p==1],max(ii)),max(ii[p==0],min(ii))))/2}
identMarkers=scan("CellIdentityMarkers.txt","")

args=commandArgs(trailing=T)
doExclusions=TRUE
i=1
dd=loadHaloObjFile(args[i],exclude=doExclusions)

dx=dd %>%
    select(Sample,SPOT,UUID,Marker,ValueType,Value) %>%
    spread(ValueType,Value)

thresholdTbl=dx %>%
    group_by(Sample,Marker) %>%
    summarize(Theta=thresholdIntensity(Intensity,Positive)) %>%
    unite(SampleMarker,Sample,Marker,sep=":",remove=F)

thetas=thresholdTbl$Theta
names(thetas)=thresholdTbl$SampleMarker

dx$Threshold=thetas[paste(dx$Sample,dx$Marker,":")]

intensityStats=dx %>%
    group_by(Sample,Marker) %>%
    summarize(Median=median(Intensity),
        Mean=mean(Intensity),
        Max=max(Intensity),
        MaxN=max(Intensity*(1-Positive)),
        MinP=min(ifelse(Positive,Intensity,Inf)))
intensityStats=full_join(thresholdTbl,intensityStats)

studyName="MelanomaIL2_Final"

markers=dx %>% distinct(Marker) %>% pull(Marker)

spots=dx %>% distinct(SPOT) %>% pull(SPOT)
cutSpots=spots[floor(len(spots)/2)]+.5
markers=dx %>% distinct(Marker) %>% pull(Marker)
samples=dx %>% distinct(Sample) %>% pull(Sample)

pltFile=cc("NewThreshold/thresholdROCCurvesSOX10",paste0(samples,collapse=","),".pdf")
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
            filter(Marker %in% identMarkers)

        markerPos="SOX10"

        markersToTest=tbl %>% filter(!Marker %in% c("SOX10","CD3")) %>% top_n(3) %>% pull(Marker)

        markersToTest=c(markersToTest,"CD3")

        #par(mfrow=c(2,3))
        layout(rbind(c(1,2,3),c(4,5,5)))
        par(pty='s')

        for(markerNeg in markersToTest) {
            ss=try({d.roc=getROC(dd %>% filter(Sample==sample & SPOT==spot),markerPos,markerNeg)})
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
