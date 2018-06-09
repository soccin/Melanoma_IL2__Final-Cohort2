require(tidyr)
require(tibble)
require(readr)

source("Halo/loadHaloObjectFiles.R")

thresholdIntensity<-function(ii,p){sum(c(min(ii[p==1],max(ii)),max(ii[p==0],min(ii))))/2}

args=commandArgs(trailing=T)
doExclusions=TRUE

dd=loadHaloObjFile(args[1],exclude=doExclusions)

sample=unique(dd$Sample)

markerToReset="SOX10"
unreasonableCombinations=scan("UnreasonableCombinations.txt","")
unreasonableCombinations=unreasonableCombinations[grepl(markerToReset,unreasonableCombinations)]
exclusionMarkers=setdiff(unique(unlist(flatten(strsplit(unreasonableCombinations,",")))),markerToReset)

dCell=dd %>%
    filter(Marker %in% c(markerToReset,exclusionMarkers)) %>%
    select(-Marker,-ValueType,-Value,-mUUID,-Object_Id) %>%
    distinct

dx=dd %>%
    filter(Marker %in% c(markerToReset,exclusionMarkers)) %>%
    select(Sample,UUID,Marker,ValueType,Value)

thresholdTbl=dx %>%
    spread(ValueType,Value) %>%
    group_by(Sample,Marker) %>%
    summarize(Theta=thresholdIntensity(Intensity,Positive)) %>%
    unite(SampleMarker,Sample,Marker,sep=":",remove=F)

thetas=thresholdTbl$Theta
names(thetas)=thresholdTbl$SampleMarker

dx=dx %>%
    select(-Sample) %>%
    unite(MarkerValueType,Marker,ValueType) %>%
    spread(MarkerValueType,Value)

dx[[cc(markerToReset,"Positive.Orig")]]=dx[[cc(markerToReset,"Positive")]]

validUUIDs=dx %>%
    select(-matches(markerToReset)) %>%
    select(matches("Positive|UUID")) %>%
    mutate(AnyPos=rowSums(.[-1])>0) %>%
    filter(!AnyPos) %>% pull(UUID)

dx=dx %>% filter(UUID %in% validUUIDs) %>% select(matches(markerToReset),UUID)

dx=dCell %>% filter(UUID %in% validUUIDs) %>% select(UUID,Sample,SPOT) %>% full_join(dx)

newThresholdTbl=read_csv(cc("NewThreshold/rocStats",sample,"SOX10_.csv")) %>%
    filter(MarkerNeg=="CD3" & auc>.6 & precOpt>.85 & thetaOpt<thetaOrig) %>%
    select(Sample,Spot,thetaOpt) %>%
    rename(SPOT=Spot)

dReset=left_join(dx,newThresholdTbl) %>%
    filter(!is.na(thetaOpt)) %>%
    mutate(SOX10_Positive=case_when(SOX10_Intensity>thetaOpt ~ 1, T ~ SOX10_Positive.Orig)) %>%
    select(UUID,SOX10_Positive) %>%
    rename(SOX10_Positive.New=SOX10_Positive) %>%
    mutate(Marker="SOX10",ValueType="Positive")

d.Orig=dd %>%
    filter(UUID %in% dReset$UUID & Marker==markerToReset & ValueType=="Positive") %>%
    mutate(ValueType="Positive.Orig")

dd=left_join(dd,dReset) %>%
    mutate(Value=ifelse(!is.na(SOX10_Positive.New),SOX10_Positive.New,Value)) %>%
    select(-SOX10_Positive.New)

dd=bind_rows(dd,d.Orig)
newRDA=file.path("NewThreshold",gsub(".rda","__ReThreshold_SOX10.rda",basename(args[1])))
saveRDS(dd,file=newRDA,compress=T)

