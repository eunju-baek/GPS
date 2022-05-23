getmode <- function(v) {
 uniqv <- unique(v)
 uniqv[which.max(tabulate(match(v, uniqv)))]
}

G2_mode <- big_apply(G2_ht,a.FUN= function(X,ind){getmode(X[,df_beta[["_NUM_ID_"]]])},a.combine= "rbind")
