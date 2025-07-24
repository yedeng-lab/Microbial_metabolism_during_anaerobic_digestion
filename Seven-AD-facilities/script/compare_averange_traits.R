#compare averange traits of molecular composition
#between input and reactors
#@yang
rm(list=ls())

##load data
original_data=read.csv("molecular_list.csv")
group_data=read.csv("sample_list.csv")
names(group_data)[1]="id"
purple_data=original_data[group_data$id]
use_group=group_data[3]
ave=read.csv("average_traits.csv",row.names = 1)
com_ave=data.frame(matrix(ncol=8,nrow=19))
##mean1 input mean2 output
colnames(com_ave)=c("parameter","mean","sd","se","mean","sd","se","p")

for (ii in 2:20) {
  jj=ii-1
  the_data=ave[ii]
  com_ave[jj,1]=colnames(the_data)
  com_ave[jj,2]=aggregate(x=the_data[1],by=list(use_group$type),mean)[1,2]
  com_ave[jj,5]=aggregate(x=the_data[1],by=list(use_group$type),mean)[2,2]
  com_ave[jj,3]=aggregate(x=the_data[1],by=list(use_group$type),sd)[1,2]
  com_ave[jj,6]=aggregate(x=the_data[1],by=list(use_group$type),sd)[2,2]
  com_ave[jj,4]=aggregate(x=the_data[1],by=list(use_group$type),FUN = function(x) sd(x)/sqrt(length(x)))[1,2]
  com_ave[jj,7]=aggregate(x=the_data[1],by=list(use_group$type),FUN = function(x) sd(x)/sqrt(length(x)))[2,2]
  test_result=pairwise.wilcox.test((the_data[1][,]),use_group$type,p.adjust.method = "none")
  com_ave[jj,8]=test_result$p.value
}
library(tidyverse)
marker=com_ave[8] %>% mutate(marker = cut(p, breaks = c(0, 0.001, 0.01, 0.05,Inf),
                                          labels = c("***", "**", "*","")))
com_ave$marker=marker$marker

