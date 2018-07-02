reassignAndPlotSampleSpot<-function(dd) {

    dx <- spreadMarkerTbl(dd)
    dx <- reassignSOX10(dx)

    dCell <- getCellTable(dd)
    dCell <- full_join(dCell,dx)

    plotSpatialCellReassign(dCell)

}

source("toolsReAssignDiag.R")
source("Halo/loadHaloObjectFiles.R")

targetMarker="SOX10"
positiveMarkers=c("S100B")
identMarkers=getIdentityMarkers()
negativeMarkers=setdiff(identMarkers,c(targetMarker,positiveMarkers))
allMarkers="CD20:CD8:CD3:PCK26"
thresholdDir="NewThresholdV5r1/SOX10_S100B/CD20,CD8,CD3,PCK26"

stats <- dir_ls("NewThresholdV5r1",glob="*.csv") %>%
    map(read_csv) %>%
    bind_rows %>%
    select(-PCT) %>%
    spread(State,n) %>%
    mutate(PCT.Delta=100*ifelse(is.na(`1>0`),`0>1`/`1>1`,`1>0`/`1>1`)) %>%
    mutate(DAPI=rowSums(.[3:6],na.rm=T)) %>%
    mutate(PCT.SOX10.Orig=100*`1>1`/DAPI)


lowSox10=stats %>% arrange(desc(PCT.Delta)) %>% filter(PCT.SOX10.Orig<10 & PCT.Delta>10)
fovsToPlot=bind_rows(lowSox10)
highDeltaPos=stats %>% arrange(desc(PCT.Delta)) %>% filter(PCT.SOX10.Orig>=10 & PCT.Delta>50)
fovsToPlot=bind_rows(highDeltaPos,fovsToPlot)
lowDeltaPos=stats %>% arrange((PCT.Delta)) %>% filter(PCT.Delta< -2)
fovsToPlot=bind_rows(lowDeltaPos,fovsToPlot)
fovsToPlot=fovsToPlot %>% distinct()

ddCache=list()

pdf(file="fovsToCheckReThresholdV5r1.pdf",width=11,height=8.5)

for(df in split(fovsToPlot,seq_len(nrow(fovsToPlot)))) {
    sampleID=df$Sample
    spot=df$Spot
    if(is.null(ddCache[[sampleID]])) {
        #rdaFile=paste0("../SubSampleCohort2/",sampleID,"_MegaTableV5b_Excl_Small.rda")
        rdaFile=paste0("data/byrne/",sampleID,"_MegaTableV5b_Excl.rda")
        cat("Loading",rdaFile,"\n")
        ddCache[[sampleID]]=loadHaloObjFile(rdaFile,TRUE)
    }
    dd=ddCache[[sampleID]]
    newThetas=getNewThetaTable(sampleID,targetMarker,allMarkers)
    reassignAndPlotSampleSpot(dd %>% filter(SPOT==spot))

}

dev.off()


