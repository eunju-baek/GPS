library(bigstatsr)
library(bigsnpr)
library(dplyr)
library(data.table)
library(magrittr)
library(R.utils)


obj_vali <- snp_attach("/BiO/enju07/4_LDpred/2_GPS/phase3/1_validation/mer_vali_v3.rds")
G_vali <- obj_vali$genotypes
CHR_vali <- obj_vali$map$chromosome
POS_vali <- obj_vali$map$physical.pos
G2_vali <- snp_fastImputeSimple(G+vali)



# match SNPs between sumstats and validation set
map_vali <- obj_vali$map[-(2:3)]
names(map_vali) <- c("chr", "pos", "a0", "a1")
info_snp_vali <- snp_match(sumstats2, map_vali)


tmp <- tempfile(tmpdir = "/BiO/enju07/4_LDpred/2_GPS/phase3/1_validation/tmp-data")

fam.order <- as.data.table(obj_vali$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))

phenotype <- fread("/BiO/enju07/4_LDpred/2_GPS/phase3/phenopheno_275809_vali.txt")

pheno_vali <- phenotype[fam.order, on = c("FID", "IID")]


pred_grid_vali <- big_prodMat( G2_vali, beta_grid, ind.col = info_snp$`_NUM_ID_`)



# or you can refer following address to reduce computing time[https://privefl.github.io/bigsnpr/articles/LDpred2.html]

# best_beta_grid <- grid.param %>%
#  mutate(id = row_number()) %>%
  # filter(sparse) %>% 
#  arrange(desc(score)) %>%
#  slice(1) %>%
#  pull(id) %>%
#  beta_grid[, .]


pred_vali <- big_prodVec(G, best_beta_grid, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
