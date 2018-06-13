mm=scan("reThresholdMarkers__SOX10_.txt.orig","")
mmComb=c()
for(ii in seq(len(mm))) {

    mmComb=c(mmComb,apply(combn(mm,ii),2,function(x){paste(x,collapse=":")}))

}
write(mmComb,"allCombinations")


