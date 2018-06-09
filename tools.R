getCellTypeTable<-function(){
    cellTypes=read_xlsx("meta/MelanomaCohort2_CellTypes.xlsx")
    cellTypes$MarkersNormalized=strsplit(cellTypes$Marker_combination,",") %>%
        map(function(x){paste0(sort(x),collapse=",")}) %>%
        unlist
    cellTypes
}
