
spreadMarkerTbl <- function(din) {
    din %>%
        filter(Marker %in% identMarkers) %>%
        select(Sample,SPOT,UUID,Marker,ValueType,Value) %>%
        rename(Spot=SPOT) %>%
        unite(MVT,Marker,ValueType) %>%
        spread(MVT,Value) %>%
        select(-matches("_Intensity")) %>%
        select(cc(negativeMarkers,"Positive")) %>%
        apply(.,1,function(x){all(x==0)}) %>%
        select(-matches("_Positive"),SOX10_Positive,superNeg)
}

getNewThetaTable <- function(sampleID,targetMarker,allNegMarkers) {
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

reassignSOX10 <- function(dx,RULE=1) {
    if(RULE==1) {
        ODIR="NewThresholdV5"
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
    dx %>% mutate(State=paste0(SOX10_Positive.Orig,">",SOX10_Positive))
}

getCellTable <- function(dd) {
    dd %>%
        select(UUID,Sample,XMin,XMax,YMin,YMax,SPOT) %>%
        distinct %>%
        mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2)
}

plotSpatialCellReassign <- function(dCell) {

        PADDING=30
        CLIPFOV=TRUE
        clipSize=PADDING
        fovTag=paste0("clip",PADDING)
        bbData=padBoundingBox(bbFOV0,-clipSize/pixel2um)

        bbPlot=bbFOV0
        bbPlot$X0=bbPlot$X0
        bbPlot$Y1=bbPlot$Y1

        colThres=c("#BEBEBE","#BEBEBE",
                "#17d479","#1772d4",
                "#d41a17","#d41a17",
                "#e3a19e","#e3a19e")

        cexThres=c(.5,.5,1,1,1,1,.75,.75)

        sampleName <- dCell %>% distinct(Sample) %>% pull
        spot <- dCell %>% distinct(SPOT) %>% pull

        ds=dCell %>%
            mutate(reThresFlag=2*SOX10_Positive.Orig+4*SOX10_Positive+ifelse(superNeg,1,0)+1)

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
        drawBoundariesForSample(sampleName)

        #########################################################################################
        #
        # Customize cell type/marker points
        i.pos=(ds$reThresFlag>1)

        with(ds[i.pos,],
            points(X,Y,pch=16,col=colThres[reThresFlag],cex=cexThres[reThresFlag])
            )

        #########################################################################################

}

# newThetas=getNewThetaTable(sampleID,targetMarker,allMarkers)
# dx <- spreadMarkerTbl(dd)
# dx <- reassignSOX10(dx)

# dCell <- getCellTable(dd)
# dCell <- full_join(dCell,dx)
