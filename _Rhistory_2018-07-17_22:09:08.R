require(tidyverse)
require(readxl)
dd=read_xlsx("cellTypeCombinationTable__reThresV5r2.xlsx")
dd
d1=dd
d1$Markers
gsub("CD68","",d1$Markers)
gsub("CD68","",d1$Markers) %>% gsub("^,|,,|,$","",.)
gsub("CD68","",d1$Markers) %>% gsub("^,|,,|,$","",.) %>% len
gsub("CD68","",d1$Markers) %>% gsub("^,|,,|,$","",.) %>% uniq(.) %>% len
gsub("CD68","",d1$Markers) %>% gsub("^,|,,|,$","",.) %>% unique(.) %>% len
gsub("CD68","",d1$Markers) %>% gsub("^,|,,|,$","",.)
d1$Markers=gsub("CD68","",d1$Markers) %>% gsub("^,|,,|,$","",.)
d1
d1 %>% group_by(Markers)
d1 %>% group_by(Markers) %>% summarise_if(is.numeric,sum)
d1 %>% group_by(Markers) %>% summarise_if(is.numeric,sum,na.rm=T)
d1.collapse=d1 %>% group_by(Markers) %>% summarise_if(is.numeric,sum,na.rm=T)
dd %>% filter(grepl("CD68",Markers))
dd %>% filter(!grepl("CD68",Markers))
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers")
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total"))
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|mel_2"))
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|mel_2."))
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|mel_2\\."))
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|mel_2\\.|Marker"))
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker"))
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker")) %>% mutate(Error=(Total.y-Total.x)/Total.x)
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker")) %>% mutate(Error=100*(Total.y-Total.x)/Total.x)
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker")) %>% mutate(PCT.Error=100*(Total.y-Total.x)/Total.x)
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker")) %>% mutate(PCT.Error=100*(Total.y-Total.x)/Total.x) %>% arrange(desc(PCT.Error)
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker")) %>% mutate(PCT.Error=100*(Total.y-Total.x)/Total.x) %>% arrange(desc(PCT.Error))
dd
dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker|Cell_")) %>% mutate(PCT.Error=100*(Total.y-Total.x)/Total.x) %>% arrange(desc(PCT.Error))
error=dd %>% filter(!grepl("CD68",Markers)) %>% left_join(d1.collapse,by="Markers") %>% select(matches("Total|Marker|Cell_")) %>% mutate(PCT.Error=100*(Total.y-Total.x)/Total.x) %>% arrange(desc(PCT.Error))
require(openxlsx)
write.xlsx(error,"pctErrorForMissing_CD68.xlsx")
