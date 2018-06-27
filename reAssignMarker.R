require(tidyverse)
require(fs)

source("getCellTypes.R")
source("Halo/loadHaloObjectFiles.R")

args <- commandArgs(trailing=T)
dd <- loadHaloObjFile(args[1],exclude=TRUE)
sampleID=dd %>% distinct(Sample) %>% pull

if(len(sampleID)!=1) {
    cat("can only process 1 sample\n")
    quit()
}

targetMarker="SOX10"
positiveMarkers=c("S100B")
identMarkers=getIdentityMarkers()
negativeMarkers=setdiff(identMarkers,c(targetMarker,positiveMarkers))
allMarkers="CD20:CD8:CD3:PCK26"

getNewThetaTable <- function(sampleID,targetMarker,allNegMarkers) {
    thresholdDir="NewThresholdV4/SOX10_S100B/CD20,CD8,CD3,PCK26"
    thetaFile=file.path(thresholdDir,cc("rocStats",sampleID,targetMarker,".csv"))
    xx=read_csv(thetaFile)
    baseTbl=xx %>% select(Sample,Spot,thetaOrig,numDAPI,numPOS) %>% distinct()
    reSetThres=xx %>%
        filter(MarkerNeg==allMarkers & auc>.8 & precOpt>.65 & numPOS/numDAPI>.02) %>%
        select(Sample,Spot,thetaOpt,auc,numDAPI,numPOS,PCT.Delta)
    newThetaTable=left_join(baseTbl,reSetThres) %>%
        mutate(thetaNew=ifelse(is.na(thetaOpt),thetaOrig,thetaOpt))
    newThetaTable
}

newThetas=getNewThetaTable(sampleID,targetMarker,allMarkers)

dx=dd %>%
    filter(Marker %in% identMarkers) %>%
    select(Sample,SPOT,UUID,Marker,ValueType,Value) %>%
    rename(Spot=SPOT) %>%
    unite(MVT,Marker,ValueType) %>%
    spread(MVT,Value) %>%
    select(-matches("_Intensity"),SOX10_Intensity)

dx$superNeg=dx %>% select(cc(negativeMarkers,"Positive")) %>% apply(.,1,function(x){all(x==0)})
dx=dx %>% select(-matches("_Positive"),SOX10_Positive)

dx <- left_join(dx,newThetas %>%
    select(Sample,Spot,thetaNew)) %>%
    mutate(SOX10_PositiveNew=case_when(
        SOX10_Intensity>thetaNew & superNeg ~ 1,
        SOX10_Intensity>thetaNew ~ SOX10_Positive,
        T ~ 0))

tbl=dx %>% count(SOX10_Positive,SOX10_PositiveNew,superNeg,SOX10_Intensity>thetaNew)

stop("BREAK")

dCell=dd %>%
    filter(Marker %in% identMarkers) %>%
    select(-Marker,-ValueType,-Value,-mUUID,-Object_Id) %>%
    distinct





dirs=dir_ls("NewThresholdV4",recursive=T,regexp="roc.*csv")

stats <- dirs[!grepl("CD45",dirs)] %>%
    map(read_csv) %>%
    bind_rows %>%
    group_by(Sample,Spot) %>%
    mutate(Rank=rank(-auc)) %>%
    ungroup

nCases=stats %>% distinct(dblPosMarker,MarkerNeg) %>% nrow(.)


stats %>% group_by(Sample,Spot) %>% mutate(Rank=rank(-auc)) %>% select(-matches("Theta|Opt")) %>% arrange(PCT.Delta) %>% filter(Sample=="mel_5" & Spot==2) %>% arrange(Rank) %>% data.frame
