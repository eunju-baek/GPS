setwd(getwd())

library(bigstatsr)
library(bigsnpr)
library(dplyr)
library(data.table)
library(magrittr)
library(R.utils)

## Load the Rdata containing the correlation matrix calculated in Step 2 (quantitative trait GPS)
load("Score.Rdata")

sum_t2d <- fread("/home/adam/BiO/Hyein/90Traits/T2D/input_GPS_T2D_v1.txt")
> str(sum_t2d)
Classes ??data.table?? and 'data.frame':	1345538 obs. of  9 variables:
 $ chr    : int  1 1 1 1 1 1 1 1 1 1 ...
 $ rsid   : chr  "rs12565286" "rs3094315" "rs3131972" "rs3115860" ...
 $ pos    : int  721290 752566 752721 753405 754182 760912 761147 761732 768448 776546 ...
 $ a0     : chr  "C" "G" "A" "C" ...
 $ a1     : chr  "G" "A" "G" "A" ...
 $ beta   : num  -0.0402 0.0189 0.0178 0.0244 0.0233 ...
 $ beta_se: num  0.0444 0.0214 0.0213 0.0229 0.0229 ...
 $ p      : num  0.365 0.376 0.403 0.287 0.307 ...
 $ n_eff  : int  34881 34881 34881 34881 34881 34881 34881 34881 34881 34881 ...
 - attr(*, ".internal.selfref")=<externalptr> 
 
> info_snp_t2d <- snp_match(sum_t2d,map_ht)
1,345,538 variants to be matched.
0 ambiguous SNPs have been removed.
1,141,242 variants have been matched; 0 were flipped and 922 were reversed.

> df_beta_t2d <- info_snp_t2d[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]

> ldsc_t2d <- snp_ldsc( ld, length(ld), chi2 = (df_beta_t2d$beta / df_beta_t2d$beta_se)^2, sample_size = df_beta_t2d$n_eff, blocks = NULL)

> h2_t2d <- ldsc_t2d[["h2"]]

> beta_inf_t2d <- snp_ldpred2_inf(corr,df_beta_t2d,h2=h2_t2d)
pred_inf_t2d <- big_prodVec(G2_ht, beta_inf_t2d, ind.row = ind.test, ind.col = info_snp_t2d$`_NUM_ID_`)

