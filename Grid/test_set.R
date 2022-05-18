# test set for selecting grid hyper-parameter

library(bigstatsr)
library(bigsnpr)
library(dplyr)
library(data.table)
library(magrittr)
library(R.utils)

snp_readBed("/BiO/enju07/4_LDpred/2_GPS/phase3/0_test/v3/mer_test_v3.bed")
obj.bigSNP <- snp_attach("/BiO/enju07/4_LDpred/2_GPS/phase3/0_test/v3/mer_test_v3.rds")


G_obj <- obj.bigSNP$genotypes
CHR_obj <- obj.bigSNP$map$chromosome
POS_obj <- obj.bigSNP$map$physical.pos
G2_obj <- snp_fastImputeSimple(G_obj)

sumstats <- bigreadr::fread2("/BiO/jionekang/GIANT,UKB_meta_gwas_bmi/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt")
sumstats2 <- sumstats [ , c(1,3,2,4,5,7,8,9,10)]
names(sumstats2) <- c("chr","rsid","pos","a0","a1","beta","beta_se","p","n_eff")

map_obj <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp_obj <- snp_match(sumstats2, map)

phenotype <- fread("/BiO/enju07/4_LDpred/2_GPS/phase3/pheno/pheno_68952_test.txt")

fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc( ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]

y <- phenotype[fam.order, on = c("FID", "IID")]


p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2) 
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4) 
grid.param <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)) 

beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
pred_grid <- big_prodMat( G2, beta_grid, ind.col = info_snp$`_NUM_ID_`)


ind.val <- phenotype[,1]


## Hyper-parameter selection: Best R-squared

grid.param$score <- apply(pred_grid[ind.val,],2,function(x){
+  summary(lm(phenotype[ind.val,3] ~ x))$adj.r.squared
+ })



  
library(dplyr)
grid.param %>%
  mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
  arrange(desc(score)) %>%
  mutate_at(c("score", "sparsity"), round, digits = 3) %>%
  slice(1:10)  

##       p     h2 sparse score sparsity  id
## 1  0.180 0.2054  FALSE 0.308        0  48
## 2  0.320 0.2054  FALSE 0.307        0  49
## 3  0.100 0.2054  FALSE 0.307        0  47
## 4  0.560 0.2054  FALSE 0.306        0  50
## 5  1.000 0.2054  FALSE 0.305        0  51
## 6  1.000 0.2054   TRUE 0.305        0 102
## 7  0.056 0.2054  FALSE 0.301        0  46
## 8  0.100 0.1467  FALSE 0.293        0  30
## 9  0.180 0.1467  FALSE 0.292        0  31
## 10 0.320 0.1467  FALSE 0.291        0  32


## Plotting for hyper-parameter selection
library(ggplot2)
ggplot(grid.param, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "adj R squared", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
