getROC <- function(df,markerPos,markerNeg) {

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

getROCMulti <- function(df,markerPos,markerNegs) {

    dt=df %>%
        filter(Marker %in% c(markerPos,markerNegs)) %>%
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

    d.roc=roc(dt$category,dt$posIntensity)
}

plot.roc.1 <- function(d.roc,markerPos,markerNeg) {

    plot(d.roc,main=paste0(markerPos,", ",markerNeg),asp='s',
        print.thres=T,print.thres.cex=1.12,
        print.auc=T,print.auc.cex=1.12)

}

thresholdIntensity<-function(ii,p){sum(c(min(ii[p==1],max(ii)),max(ii[p==0],min(ii))))/2}

