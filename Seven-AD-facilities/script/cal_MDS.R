#analysis of Molecular dissimilarity in segments (MDS)
#@yang 20250330

rm(list=ls())
library(vegan)

##load data
original_data=read.csv("molecular_list.csv")
group_data=read.csv("sample_list.csv")
names(group_data)[1]="id"

data=original_data[c("molecular_weight",group_data$id[group_data$type=="Output"])]

##function
step_dis=function (data1,from1=100,to1=800,step1=1,halfweidth=20,method="bray",binary=T){
  loopmain=seq((from1+halfweidth),(to1-halfweidth),step1)
  loopsave=data.frame(loopmain,dissimilarity=0)#save result
  len=length(loopmain)
  for (ii in 1:len){
    para1=loopmain[ii]-halfweidth
    para2=loopmain[ii]+halfweidth
    dothework=data1[data1[1]>=para1&data1[1]<=para2,2:3]#select data
    workresult=as.matrix(vegdist(t(dothework),method =method,binary = binary))#in this method D=1-S,cal dissimilarity
    loopsave[ii,2]=workresult[1,2] 
  }
  return(loopsave)
}


##calculate MDS
jj=2
kk=2
data1=data[c(1,jj,kk)]
all_reault=step_dis(data1=data1,from1=100,to1=800,step1=5,halfweidth=20,method="bray",binary=F)
colnames(all_reault)[2]=paste(colnames(data)[jj],colnames(data)[kk],sep="_")

for (jj in 2:42){
  for (kk in (jj+1):43) {
    data1=data[c(1,jj,kk)]
    one_reault=step_dis(data1=data1,from1=100,to1=800,step1=5,halfweidth=20,method="bray",binary=F)
    colnames(one_reault)[2]=paste(colnames(data)[jj],colnames(data)[kk],sep="_")
    all_reault=cbind(all_reault,one_reault[2]) 
  }
}
all_reault2=all_reault[-2]
all_dis_weighted=all_reault2

all_dis_weighted_melt=reshape2::melt(all_dis_weighted,id.vars="loopmain")

cor.test(all_dis_weighted_melt$loopmain,all_dis_weighted_melt$value,method="spearman")

all_a=all_dis_weighted[1:2]
names(all_a)[2]="average value"
all_a$`average value`=rowMeans(all_dis_weighted[2:ncol(all_dis_weighted)],na.rm = T)


##detect change point
library(segmented)

fit_lm <- lm(average.value ~loopmain, data = step_dissimilarity)#aic -81.6471
model1=segmented(fit_lm, average.value = ~loopmain, npsi = 1)#aic -320.193
model2=segmented(fit_lm, average.value = ~loopmain, npsi = 2)#aic -347.261
model3=segmented(fit_lm, average.value = ~loopmain, npsi = 3)#aic -636.448
model4=segmented(fit_lm, average.value = ~loopmain, npsi = 4)#aic -643.991

AIC(fit_lm)
AIC(model1)
AIC(model2)
AIC(model3)
AIC(model4)

plot(model3)
