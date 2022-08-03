library(data.table)
library(dplyr)

plc <- fread("IS_stroke_UKB_WBU.csv.gz")
plc <- data.frame(plc)
hg <- fread("../Asthma/hg19_avsnp147.txt.gz")
hg <- data.frame(hg)

plc1 <- plc[,c(2,3,4,8,5,6)]


hg['position'] <- paste(hg$chrom,hg$pos_st,sep="_")
plc1['position'] <- paste(plc1$chrom,plc1$pos,sep='_')
plc2 <- left_join(plc1,hg,by='position')


plc3 <- plc2[,c(1,14,2,3,4,5,6,7)]
names(plc3) <- c("CHR","SNP","BP","A1","A2","BETA","SE","P")

fwrite(plc3,'input_IS_stroke_UKB_WBU.txt',quote=F,sep='\t',row.names=F,na="NA")
