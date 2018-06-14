# getROC <- function(df,markerPos,markerNeg) {

#     dt=df %>%
#         filter(Marker %in% c(markerPos,markerNeg)) %>%
#         unite(MarkerValue,Marker,ValueType) %>%
#         select(Sample,SPOT,UUID,MarkerValue,Value) %>%
#         spread(MarkerValue,Value)

#     posIntensity=cc(markerPos,"Intensity")
#     posPositive=cc(markerPos,"Positive")
#     negIntensity=cc(markerNeg,"Intensity")
#     negPositive=cc(markerNeg,"Positive")

#     dt[["posIntensity"]]=dt[[posIntensity]]
#     dt[["posPositive"]]=dt[[posPositive]]
#     dt[["negIntensity"]]=dt[[negIntensity]]
#     dt[["negPositive"]]=dt[[negPositive]]

#     dt %<>%
#         mutate(category=case_when(
#                 negPositive==1 ~ "NEG",
#                 TRUE ~ "POS"))

#     d.roc=roc(dt$category,dt$posIntensity)
# }

getROCMulti <- function(df,markerPos,markerNegsStr,dblPosMarker=NULL) {

    cat("markerNegsStr =",markerNegsStr,"\n")
    markerNegs=strsplit(markerNegsStr,":")[[1]]
    dt=df %>%
        filter(Marker %in% c(markerPos,markerNegs,dblPosMarker)) %>%
        unite(MarkerValue,Marker,ValueType) %>%
        select(Sample,SPOT,UUID,MarkerValue,Value) %>%
        spread(MarkerValue,Value)

    posIntensity=cc(markerPos,"Intensity")
    posPositive=cc(markerPos,"Positive")
    negIntensity=cc(markerNegs,"Intensity")
    negPositive=cc(markerNegs,"Positive")

    dt[["posIntensity"]]=dt[[posIntensity]]
    dt[["posPositive"]]=dt[[posPositive]]
    dt[["negIntensity"]]=apply(dt[negIntensity],1,max)
    dt[["negPositive"]]=ifelse(apply(dt[negPositive]==1,1,any),1,0)

    dt %<>%
        mutate(category=case_when(
                negPositive==1 ~ "NEG",
                TRUE ~ "POS"))

    if(!is.null(dblPosMarker)) {
        singlePosIdx=dt$category=="POS" & dt[[cc(dblPosMarker,"Positive")]]==0
        dt$category[singlePosIdx] <- "NEG"
    }

    d.roc=roc(dt$category,dt$posIntensity)

}

plot.roc.1 <- function(d.roc,markerPos,markerNeg) {

    plot(d.roc,main=paste0(markerPos,", ",markerNeg),asp='s',
        print.thres=T,print.thres.cex=1.12,
        print.auc=T,print.auc.cex=1.12)

}

thresholdIntensity <- function(ii,p){sum(c(min(ii[p==1],max(ii)),max(ii[p==0],min(ii))))/2}

getROCStats <- function(sample,spot,markerPos,markerNegs,d.roc) {

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
        acc0=mean(coords(d.roc,theta0,"t","acc")),
        numCases=len(d.roc$cases),
        numControls=len(d.roc$controls)
    )

    stats$F1=2*stats$precOpt*stats$sensOpt/(stats$precOpt+stats$sensOpt)
    stats$F0=2*stats$prec0*stats$sens0/(stats$prec0+stats$sens0)

    return(stats)

}