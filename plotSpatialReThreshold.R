
# suppressPackageStartupMessages(require(RColorBrewer))
# suppressPackageStartupMessages(require(gplots))
# suppressPackageStartupMessages(require(shades))

# source("Halo/funcs.R")

data(melanomaNewV20_CellInfo)

source("HaloSpatial/spatialV3.R")
source("HaloSpatial/plot.R")
source("HaloSpatial/halo.R")

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

dCell=filterCellsWithinBoundryLayer(dCell,PADDING)

annote=read_csv("data/melanomaNewV20_FOVAnnotations.csv")
fovsToRedact=annote %>% filter(EXCLUDE) %>% unite("SampFOV",Sample,FOV) %>% pull(SampFOV)

dCell=redactFOVs(dCell,fovsToRedact)

#########################################################################################
#########################################################################################
# Join cell level data to dCell for plotting
data(melanomaNewV20__ReThreshold_SOX10)

dCell=dd %>%
    select(UUID,Marker,ValueType,Value) %>%
    filter(Marker=="SOX10" & grepl("Positive",ValueType)) %>%
    unite(MarkerValueType,Marker,ValueType) %>%
    spread(MarkerValueType,Value) %>%
    right_join(dCell)

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

        plotFOVRaw(bbData,sampleName,spot,bbPlot)
        points(ds$X,ds$Y,pch=16,col=8,cex=.5)
        drawBoundariesForSample(sampleName)

        #########################################################################################
        #
        # Customize cell type/marker points
        i.pos=(ds$SOX10_Positive==1)
        i.newPos=which(ds$SOX10_Positive.Orig==0 & ds$SOX10_Positive==1)

        points(ds$X[i.pos],ds$Y[i.pos],pch=16,col="#FFC7C6",cex=.5)
        points(ds$X[i.newPos],ds$Y[i.newPos],pch=16,col="#FF0000",cex=1)
        #########################################################################################


    }

    if(!interactive()) dev.off()

}
