library(olsrr)
library(data.table)
library(dplyr)

a <- fread("pheno.txt")

a <- data.frame(a)

for(i in 1:5){
assign(paste("qq",i,sep=""),a[a$g48_5==i,])
}



 for(i in 1:5){assign(paste("mod",i,sep=''), paste('qq',i,sep='') %>% get() %>% lm(bmi ~ V48 + + X1558.0.0 + (X1558.0.0*V48)+age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+array,data=.)%>% ols_test_score())}
 
 mat <- matrix(nrow=5,ncol=2)
 
 for(i in 1:5){mat[i,1] <- get(paste("mod",i,sep=''))$score}
 for(i in 1:5){mat[i,2] <- get(paste("mod",i,sep=''))$p}
 
