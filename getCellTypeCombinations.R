source("Halo/loadHaloObjectFiles.R")
source("getCellTypes.R")

args=commandArgs(trailing=T)
markers=getIdentityMarkers()

doExclusions=TRUE
i=1
dd=loadHaloObjFile(args[i],exclude=doExclusions)

tbl=dd %>%
    filter(ValueType=="Positive" & Marker %in% markers) %>%
    select(UUID,Sample,Marker,Value) %>%
    spread(Marker,Value) %>%
    #mutate_all(funs(replace(., is.na(.), 0))) %>%
    select(-UUID) %>%
    group_by_at(c("Sample",markers)) %>%
    summarise(Count=n()) %>%
    arrange(desc(Count)) %>%
    ungroup

tbl$MarkerSig=tbl %>%
    select(-Sample,-Count) %>%
    apply(.,1,function(x){paste0(markers[x==1],collapse=",")})

#num.DAPI=dd %>% filter(Marker=="DAPI" & ValueType=="Positive") %>% nrow

sample=tbl$Sample[1]

if(doExclusions) {
    write_csv(tbl,cc("cellTypeCombinations",sample,"ReThreshold.csv"))
} else {
    stop("STOP DOING THIS")
    #write_csv(tbl,cc("cellTypeCombinations","NoExclusions",sample,".csv"))
}

