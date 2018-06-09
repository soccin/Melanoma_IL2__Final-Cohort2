require(tidyverse)
require(readxl)
require(fs)

getIdentityMarkers<-function() {
    markerMetaFile=dir_ls("meta",regex="/[A-Za-z].*_MarkerDescriptions.xlsx")
    read_xlsx(markerMetaFile) %>%
        filter(Cell_Type=="X") %>%
        pull(Marker_name) %>%
        sort
}

getCellTypes<-function(dd) {

    cellIdentityMarkers=getIdentityMarkers()

    dm=dd %>%
        filter(ValueType=="Positive" & Marker %in% cellIdentityMarkers) %>%
        select(Sample,SPOT,UUID,Marker,Value) %>%
        mutate(Marker=paste0("mrk_",Marker)) %>%
        spread(Marker,Value) %>%
        select(matches("UUID|mrk_")) %>%
        data.frame %>%
        column_to_rownames("UUID")

    colnames(dm)=gsub("mrk_","",colnames(dm))

    cellTypeStr=unlist(
        pmap(dm,
            function(...){
                xx=list(...);
                paste0(sort(names(xx)[xx==1]),collapse=",")
            }
            )
        )

    cellTypes=as.tibble(cbind(UUID=rownames(dm),CellType=cellTypeStr))

    left_join(dd,cellTypes)

}

