require(tidyverse)
require(readxl)
require(magrittr)
require(openxlsx)

source("tools.R")

args <- commandArgs(trailing=T)

dd <- args %>%
    map(read_csv) %>%
    bind_rows %>%
    mutate(MarkerSig = ifelse(is.na(MarkerSig),"AllNeg",MarkerSig))

samples <- dd %>% distinct(Sample) %>% pull(Sample) %>% sort

dd <- dd %>%
    select(Sample,MarkerSig,Count) %>%
    spread(Sample,Count) %>%
    mutate(Total=rowSums(.[-1],na.rm=T)) %>%
    arrange(desc(Total))

cellTypes <- getCellTypeTable()

fixNames <- cellTypes$Marker_combination
names(fixNames) <- cellTypes$MarkersNormalized

cellTypeName <- cellTypes$Cell_type
names(cellTypeName) <- cellTypes$MarkersNormalized

cellSubType <- cellTypes$Subtype
names(cellSubType) <- cellTypes$MarkersNormalized

dd %<>%
    mutate(Markers=ifelse(is.na(fixNames[MarkerSig]),MarkerSig,fixNames[MarkerSig])) %>%
    mutate(Cell_Type=ifelse(is.na(cellTypeName[MarkerSig]),"",cellTypeName[MarkerSig])) %>%
    mutate(Cell_Subtype=ifelse(is.na(cellSubType[MarkerSig]),"",cellSubType[MarkerSig]))

tbl <- dd %>% select(Markers,Cell_Type,Cell_Subtype,Total,samples)

DAPI.counts <- tbl %>%
    select(samples) %>%
    summarize_all(sum,na.rm=T) %>%
    mutate(Total=rowSums(.),Markers="DAPI",Cell_Type="",Cell_Subtype="")

tbl <- bind_rows(tbl,DAPI.counts) %>%
    arrange(desc(Total))

if(grepl("NoExclusion",args[1])) {
    xFile <- "cellTypeCombinationTable_NoExclusions.xlsx"
} else {
    xFile <- "cellTypeCombinationTable.xlsx"
}

write.xlsx(tbl,xFile)
