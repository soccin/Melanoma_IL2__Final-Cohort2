require(tidyverse)
require(fs)

dirs=dir_ls("NewThresholdV4",recursive=T,regexp="roc.*csv")

stats <- dirs[!grepl("CD45",dirs)] %>%
    map(read_csv) %>%
    bind_rows %>%
    group_by(Sample,Spot) %>%
    mutate(Rank=rank(-auc)) %>%
    ungroup

nCases=stats %>% distinct(dblPosMarker,MarkerNeg) %>% nrow(.)


stats %>% group_by(Sample,Spot) %>% mutate(Rank=rank(-auc)) %>% select(-matches("Theta|Opt")) %>% arrange(PCT.Delta) %>% filter(Sample=="mel_5" & Spot==2) %>% arrange(Rank) %>% data.frame
