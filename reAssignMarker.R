source("HaloSpatial/spatialV3.R")
source("HaloSpatial/plot.R")
source("HaloSpatial/halo.R")
source("Halo/loadHaloObjectFiles.R")

source("getCellTypes.R")

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
thresholdFilesDir="NewThresholdV5r1/SOX10_S100B/CD20,CD8,CD3,PCK26"

getNewThetaTable <- function(sampleID,targetMarker,allNegMarkers) {
    thetaFile=file.path(thresholdFilesDir,cc("rocStats",sampleID,targetMarker,".csv"))
    xx=read_csv(thetaFile)
    baseTbl=xx %>% select(Sample,Spot,thetaOrig,numDAPI,numPOS) %>% distinct()
    reSetThres=xx %>%
        filter(MarkerNeg==allMarkers & auc>.8 & precOpt>.65 & numPOS/numDAPI>.02) %>%
        select(Sample,Spot,thetaOpt,auc,numDAPI,numPOS,PCT.Delta)
    newThetaTable=left_join(baseTbl,reSetThres) %>%
        mutate(thetaNew=ifelse(is.na(thetaOpt),thetaOrig,thetaOpt))
    newThetaTable
}

spreadMarkerTbl <- function(din) {
    din %>%
        filter(Marker %in% identMarkers) %>%
        select(Sample,SPOT,UUID,Marker,ValueType,Value) %>%
        rename(Spot=SPOT) %>%
        unite(MVT,Marker,ValueType) %>%
        spread(MVT,Value) %>%
        select(-matches("_Intensity"),SOX10_Intensity)
}

newThetas=getNewThetaTable(sampleID,targetMarker,allMarkers)

# dx=dd %>%
#     filter(Marker %in% identMarkers) %>%
#     select(Sample,SPOT,UUID,Marker,ValueType,Value) %>%
#     rename(Spot=SPOT) %>%
#     unite(MVT,Marker,ValueType) %>%
#     spread(MVT,Value) %>%
#     select(-matches("_Intensity"),SOX10_Intensity)

dx <- spreadMarkerTbl(dd)

dx$superNeg=dx %>% select(cc(negativeMarkers,"Positive")) %>% apply(.,1,function(x){all(x==0)})
dx=dx %>% select(-matches("_Positive"),SOX10_Positive,superNeg)

RULE=2

if(RULE==1) {
    ODIR="NewThresholdV5r1"
    dx <- left_join(dx,newThetas %>%
        select(Sample,Spot,thetaNew)) %>%
        mutate(SOX10_PositiveNew=case_when(
            SOX10_Intensity>thetaNew & superNeg ~ 1,
            SOX10_Intensity>thetaNew ~ SOX10_Positive,
            T ~ 0))
} else if(RULE==2) {
    ODIR="NewThresholdV5r2"
    dx <- left_join(dx,newThetas %>%
        select(Sample,Spot,thetaNew)) %>%
        mutate(SOX10_PositiveNew=case_when(
            SOX10_Intensity>thetaNew ~ 1,
            T ~ 0))
} else {
    stop("INVALID RULE")
}

dir.create(ODIR, showWarnings = FALSE)

tbl=dx %>% count(SOX10_Positive,SOX10_PositiveNew,superNeg,SOX10_Intensity>thetaNew)

dx %<>% rename(SOX10_Positive.Orig=SOX10_Positive,SOX10_Positive=SOX10_PositiveNew)

dgx=dx %>%
    select(Sample,Spot,UUID,matches("SOX10_")) %>%
    gather(MVT,Value,matches("SOX10_")) %>%
    separate(MVT,c("Marker","ValueType"),sep="_") %>%
    rename(SPOT=Spot,Value.New=Value) %>%
    filter(ValueType=="Positive")

dd.new=left_join(dd,dgx,by=c("UUID", "Marker", "Sample", "SPOT", "ValueType")) %>%
    mutate(Value=ifelse(!is.na(Value.New),Value.New,Value)) %>%
    select(-Value.New)

outFile=file.path(ODIR,gsub(".rda","___reThresRule1.rda",basename(args[1])))
saveRDS(dd.new,outFile,compress=T)

dx <- dx %>% mutate(State=paste0(SOX10_Positive.Orig,">",SOX10_Positive))

stats <- dx %>%
    group_by(Sample,Spot) %>%
    count(State) %>%
    mutate(n=ifelse(State=="1>0",-n,n)) %>%
    mutate(PCT=n/sum(n))

sampleName <- dx %>% distinct(Sample) %>% pull

outFile <- file.path(ODIR,cc("reThreshold_SOX10_v5",sampleName,"Stats.csv"))
write_csv(stats,outFile)

dCell = dd %>%
    select(UUID,Sample,XMin,XMax,YMin,YMax,SPOT) %>%
    distinct %>%
    mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2)

dCell=full_join(dCell,dx)

PADDING=30
CLIPFOV=TRUE
clipSize=PADDING
fovTag=paste0("clip",PADDING)
bbData=padBoundingBox(bbFOV0,-clipSize/pixel2um)

sampleNames=dCell %>% distinct(Sample) %>% pull

if(interactive()) {
    sampleNames=sampleNames[1]
}

for(sampleName in sampleNames) {

    spots=getSpotsInSample(dCell,sampleName)

    bbPlot=bbFOV0
    bbPlot$X0=bbPlot$X0
    bbPlot$Y1=bbPlot$Y1

    if(!interactive()) {
        pdf(file=file.path(ODIR,cc("reThreshold_SOX10_v5",sampleName,".pdf")),width=11,height=8.5)
    } else {
        i=1
        #stop("BREAK-A")
    }

    for(i in seq(spots)) {

        spot=spots[i]
        cat("spot=",spot,"\n")

        ds=dCell %>%
            filter(SPOT==spot & Sample==sampleName) %>%
            mutate(reThresFlag=2*SOX10_Positive.Orig+4*SOX10_Positive+ifelse(superNeg,1,0)+1)

        colThres=c("#BEBEBE","#BEBEBE",
                    "#17d479","#1772d4",
                    "#d41a17","#d41a17",
                    "#e3a19e","#e3a19e")

        cexThres=c(.5,.5,1,1,1,1,.75,.75)

        nCells=nrow(ds)
        nOrigPos=sum(ds$SOX10_Positive.Orig==1)
        nFinalPos=sum(ds$SOX10_Positive==1)
        delta=nFinalPos-nOrigPos
        rethresStats=paste0(
            "Cells=",nCells,
            " SOX10+[orig]=",nOrigPos,
            " SOX10+[rethres]=",nFinalPos,
            " delta=",round(100*delta/nOrigPos,2),"%")

        plotFOVRaw(bbData,sampleName,spot,bbPlot)
        mtext(rethresStats,3,0)

        points(ds$X,ds$Y,pch=16,col=colThres[1],cex=cexThres[1])
        drawBoundariesForSample(sampleName,spot)

        #########################################################################################
        #
        # Customize cell type/marker points
        i.pos=(ds$reThresFlag>1)

        with(ds[i.pos,],
            points(X,Y,pch=16,col=colThres[reThresFlag],cex=cexThres[reThresFlag])
            )

        #########################################################################################


    }

    if(!interactive()) dev.off()

}
