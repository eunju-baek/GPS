#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
phenotype = args[1]
environment = args[2]


library(data.table)
library(dplyr)


pheno <- fread(phenotype)

# env file is expected to be a raw("eid",...| ex. 30101..etc.) file
env <- fread(environment)

names(env)[1] <- "FID"


a <- left_join(pheno,env,by= "FID")


#you can customize filter

filt <- paste(names(env)[2:ncol(env)],collapse=">=0&")

filt <- as.logical(filt)

a <-filter(a,filt)

fwrite(a,'pheno_afterQC.txt',quote=F,sep='\t',row.names=F,na="NA")
