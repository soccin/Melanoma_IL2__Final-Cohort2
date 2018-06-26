
# suppressPackageStartupMessages(require(RColorBrewer))
# suppressPackageStartupMessages(require(gplots))
# suppressPackageStartupMessages(require(shades))

# source("Halo/funcs.R")


source("HaloSpatial/spatialV3.R")
source("HaloSpatial/plot.R")
source("HaloSpatial/halo.R")
source("Halo/loadHaloObjectFiles.R")


thresholdIntensity<-function(ii,p){sum(c(min(ii[p==1],max(ii)),max(ii[p==0],min(ii))))/2}

#
# Clean up data
# Exclude FOV's
# redact cells in padded area
#
PADDING=30
CLIPFOV=TRUE
if(CLIPFOV) {
    clipSize=PADDING
    fovTag=paste0("clip",PADDING)
    bbData=padBoundingBox(bbFOV0,-clipSize/pixel2um)
} else {
    bbData=bbFOV0
    fovTag="ORIG"
    clipSize=0
}

args <- commandArgs(trailing=T)
doExclusions <- TRUE

oDir=file.path("ReassignV4")
dir.create(oDir,showWarnings=F,recursive=T)

markerToReset="SOX10"

i=1
dd=loadHaloObjFile(args[i],exclude=doExclusions)
dd.o=dd
sampleNames=dd %>% distinct(Sample) %>% pull
if(len(sampleNames)>1) {
    cat("Can only process 1 sample at a time",len(sampleNames),"\n")
    stop("ERROR: too many samples")
}
sample=sampleNames[1]

rocStatsFile=cc("NewThresholdV4/SOX10_S100B/CD20,CD8,CD3,PCK26/rocStats",sample,markerToReset,".csv")
rocStats=read_csv(rocStatsFile) %>%
    filter(MarkerNeg=="CD20:CD8:CD3:PCK26") %>%
    filter(auc>.6 & precOpt>.65)
newThresholdTbl <- rocStats %>% select(Sample,Spot,thetaOpt) %>% rename(SPOT=Spot)

exclusionMarkers=read_xlsx("meta/MelanomaCohort2_UnreasonableCombinations_Matrix.xlsx") %>%
    select(X__1,markerToReset) %>%
    filter(SOX10=="X") %>%
    pull(X__1)

dx=dd %>%
    filter(Marker %in% c(markerToReset,exclusionMarkers)) %>%
    select(Sample,UUID,Marker,ValueType,Value)

thresholdTbl=dx %>%
    spread(ValueType,Value) %>%
    group_by(Sample,Marker) %>%
    summarize(Theta=thresholdIntensity(Intensity,Positive)) %>%
    unite(SampleMarker,Sample,Marker,sep="_",remove=F)

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

dCell = dd %>%
    select(UUID,Sample,SubSample,Layer,XMin,XMax,YMin,YMax,
        Nucleus_Area,Cytoplasm_Area,Membrane_Perimeter,Cell_Area,SPOT,SLICE) %>%
    distinct %>%
    mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2)

dCell=filterCellsWithinBoundryLayer(dCell,PADDING)

dx=dx %>% filter(UUID %in% dCell$UUID)

dCell=left_join(dCell,dx) %>% left_join(newThresholdTbl)

dCell <- dCell %>%
    mutate(SOX10_Positive=case_when(
        is.na(thetaOpt) ~ SOX10_Positive.Orig,
        !is.na(thetaOpt) & SOX10_Intensity>thetaOpt ~ 1,
        T ~ 0))

if(interactive()) stop("BREAK-A")

dx.new <- dx %>%
    select(UUID,SOX10_Positive) %>%
    gather(MVT,Value,-UUID) %>%
    separate(MVT,c("Marker","ValueType")) %>%
    rename(Value.New=Value)

dd.o <- left_join(dd.o,dx.new) %>%
    mutate(Value=case_when(!is.na(Value.New) ~ Value.New, T ~ Value)) %>%
    select(-Value.New)

newFile=file.path("ReassignV4",gsub(".rda","___reThreshV4.rda",basename(args[i])))
saveRDS(dd.o,newFile)

#########################################################################################
#########################################################################################

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
        pdf(file=cc("reThreshold_SOX10_v4",sampleName,".pdf"),width=11,height=8.5)
    } else {
        i=1
        stop("BREAK-A")
    }

    for(i in seq(spots)) {

        spot=spots[i]
        cat("spot=",spot,"\n")

        ds=dCell %>% filter(SPOT==spot & Sample==sampleName)

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

        points(ds$X,ds$Y,pch=16,col=8,cex=.5)
        drawBoundariesForSample(sampleName)

        #########################################################################################
        #
        # Customize cell type/marker points
        i.pos=(ds$SOX10_Positive==1)
        i.newPos=which(ds$SOX10_Positive.Orig==0 & ds$SOX10_Positive==1)

        points(ds$X[i.pos],ds$Y[i.pos],pch=16,col="#883333",cex=.75)
        points(ds$X[i.newPos],ds$Y[i.newPos],pch=16,col="#FF0000",cex=1)
        #########################################################################################


    }

    if(!interactive()) dev.off()

}
