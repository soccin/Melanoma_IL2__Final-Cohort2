require(tidyverse)

thresholdIntensity<-function(ii,p){sum(c(min(ii[p==1],max(ii)),max(ii[p==0],min(ii))))/2}

getThresholds<-function(dd) {
    thresholdTbl=dd %>%
        spread(ValueType,Value) %>%
        group_by(Sample,Marker) %>%
        summarize(Theta=thresholdIntensity(Intensity,Positive)) %>%
        unite(SampleMarker,Sample,Marker,sep=":",remove=F)

    thetas=thresholdTbl$Theta
    names(thetas)=thresholdTbl$SampleMarker

    return(thetas)
}
