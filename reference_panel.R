library(bigstatsr)
library(bigsnpr)
library(dplyr)
library(data.table)
library(magrittr)
library(R.utils)

# create rds file

snp_readBed("ref_mer_v3.bed") 
info <- readRDS("ref_mer_v3.rds")

G <- info$genotypes
CHR <- info$map$chromosome
POS <- info$map$physical.pos
G2 <- snp_fastImputeSimple(G)

NCORES <- nb_cores()
sumstats <- bigreadr::fread2("/BiO/jionekang/GIANT,UKB_meta_gwas_bmi/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt")
sumstats2 <- sumstats [ , c(1,3,2,4,5,7,8,9,10)]
names(sumstats2) <- c("chr","rsid","pos","a0","a1","beta","beta_se","p","n_eff")

map <- info$map[-(2:3)]


names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats2, map)

tmp <- tempfile(tmpdir = "/BiO/enju07/4_LDpred/2_GPS/phase3/reference_panel_uk10000/v3/tmp-data")
setwd('/BiO/enju07/4_LDpred/2_GPS/phase3/reference_panel_uk10000/v3/')
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL

CHR <- map$chr
chr <- map$chr
POS <- map$pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
corr0 <- NULL



for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
            G2,
            ind.col = ind.chr2,
            ncores = NCORES,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}


df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]



