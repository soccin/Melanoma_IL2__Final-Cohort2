suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(yarrr))
options( java.parameters = c("-Xss2560k", "-Xmx8g") )
suppressPackageStartupMessages(require(xlsx))

source("getCellTypes.R")
source("getThresholds.R")
source("Halo/loadHaloObjectFiles.R")

plotMarkerDensitySample <- function(dIf,marker,thresholds=NULL) {
    dIf[["MARKER__"]] <- dIf[[marker]]
    dodge <- position_dodge(width=1)
    p <- ggplot(dIf,aes(x=Sample,y=asinh(MARKER__),fill=CellType)) +
        geom_violin(position=dodge,scale="width") +
        geom_boxplot(position=dodge,width=.05,outlier.shape=NA) +
        ylab(paste0("asinh(",marker,")")) +
        ggtitle(paste(marker,"by Sample")) +
        coord_flip() +
        theme(legend.position="bottom")

    if(!is.null(thresholds)) {
        aThetas <- asinh(thresholds)
        for(ai in seq(names(aThetas)))
            p <- p+annotate('segment',alpha=0.25,size=1,
                x=ai-.495,
                xend=ai+.495,
                y=aThetas[ai],
                yend=aThetas[ai])
    }
    p
}

###############################################################################
#if(interactive()) stop("INCLUDE")
###############################################################################

args <- commandArgs(trailing=T)
cat("samples=",paste0(args,collapse=","),"\n")

dd <- args %>% map(loadHaloObjFile) %>% bind_rows

dd <- dd %>% select(Sample,SPOT,UUID,Marker,ValueType,Value)

thetas <- getThresholds(dd)

dd <- getCellTypes(dd)

dI <- dd %>%
    filter(ValueType=="Intensity") %>%
    select(-ValueType) %>%
    spread(Marker,Value) %>%
    filter(!is.na(CD25))

pMarker <- "CD25"
dI %<>% mutate(CD25.Norm=sinh(1)*CD25/thetas[paste0(Sample,":CD25")])

tCells <- scan("cellTypesTCells","")
tCells <- setdiff(tCells,"CD3,CD8,FOXP3")
mThetas <- thetas[grepl("CD25",names(thetas))]
names(mThetas) <- gsub(":.*","",names(mThetas))
norm.mThetas <- unlist(lapply(mThetas,function(x){sinh(1)}))

x11  <-  function (...) grDevices::x11(...,type='cairo')
if(interactive()) {
    x11()
} else {
#    pdf(file=paste0("markerDensity_All_",pMarker,"_.pdf"),width=11,height=8.5)
    pdir <- file.path("plots/markerDensity",gsub(".rda","",basename(args[1])))
    cat("plot directory =",pdir,"\n")
    dir_create(pdir)

    png(file=file.path(pdir,cc(pMarker,"%02d.png")),width=11,height=8.5,
        res=100,pointsize=10,type="cairo",units="in")

}


# plotMarkerDensitySample(dI %>% filter(CellType %in% tCells),"CD25",mThetas)
# plotMarkerDensitySample(dI %>% filter(CellType %in% tCells),"CD25.Norm",norm.mThetas)
# dIf %>% filter(Sample==sampleName) %>% ggplot(aes(x=factor(SPOT),asinh(CD25.Norm),fill=CellType)) + geom_violin(position=dodge,scale="width")


dIf <- dI %>% filter(CellType %in% tCells)
#yRng <- dIf %$% asinh(range(CD25.Norm))
yRng <- c(0,4.5)

# pirateplot(data=dIf,formula=asinh(CD25) ~ CellType + Sample,
#         avg.line.fun=median,inf.method="iqr",cex.names=.7,
#         main="CD25 RAW")
# legend("topleft",sort(tCells),fill=piratepal("basel",trans=.5))

# pirateplot(data=dIf,formula=asinh(CD25.Norm) ~ CellType + Sample,
#         avg.line.fun=median,inf.method="iqr",cex.names=.7,
#         main="CD25 Normalized\n[asinh(Threshold)=1]",ylim=yRng)
# legend("topleft",sort(tCells),fill=piratepal("basel",trans=.5))

samples <- dIf %>% distinct(Sample) %>% pull(Sample)

pirateplot(data=dIf,formula=asinh(CD25.Norm) ~ CellType + Sample,
        avg.line.fun=median,inf.method="iqr",cex.names=.7,
        ylim=yRng,main=paste("Samples",paste0(samples,collapse=",")))
legend("topleft",sort(tCells),fill=piratepal("basel",trans=.5))


