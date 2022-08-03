library(data.table)
library(dplyr)


fin <- fread("finngen_R6_I9_STR_EXH.gz")

hg <- fread("../Asthma/hg19_avsnp147.txt.gz")

names(hg) <- c("chrom","pos_st","pos_end","a0","a1","rsids")

fin1 <- left_join(fin,hg,by='rsids')
fin2 <- fin1[!is.na(fin1$pos_st),]
fin3 <- fin2[,c(1,5,21,4,3,9,10,7)]

names(fin3) <- c('CHR','SNP','BP','A1','A2','BETA','SE','P')



fwrite(fin3,'input_finngen_R6_I9_STR_EXH.txt',quote=F,sep='\t',row.names=F)
 
