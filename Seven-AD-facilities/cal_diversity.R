#script for calculating diversity
#@ yang
rm(list=ls())

library(vegan)

##load data
original_data=read.csv("molecular_list.csv")
group_data=read.csv("sample_list.csv")
names(group_data)[1]="id"
purple_data=original_data[group_data$id]
row.names(purple_data)=original_data$formula
group_data2=group_data[group_data$type=="Output",]



x=t(purple_data)
x=x[group_data2$id,]
Shannon <- diversity(x)
Inv_Simpson <- diversity(x, "inv")

S <- specnumber(x)
Pielou_evenness <- Shannon/log(S)
Simpson_evenness <- Inv_Simpson/S
Observed_richness <- colSums(t(x)>0)
report1 = cbind(Shannon, Inv_Simpson,Observed_richness, Pielou_evenness, Simpson_evenness)

bray.dist2=vegdist(x[group_data2$id,], method="bray",binary=F)
