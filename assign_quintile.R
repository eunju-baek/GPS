#For example, we need to assign GPS quintile group into each data frame

#root data frame >> a

a <- a %>% mutate(quintile = ntile(gps,5))

for(i in 1:5){
assign(paste0("q",i),a[a$gps==i,])
}
