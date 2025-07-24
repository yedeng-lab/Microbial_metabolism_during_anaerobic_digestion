# Variance partitioning based on multiple regression on matrices applied to segmental molecular composition differences.
#@ yang
rm(list=ls())

##load function
library(vegan)
library(ape)
library(ecodist)

### function1
###stdize from MuMIn package
source("stdize_MuMIn.R")
### function2
### DO vpa after MRM
source("RFunctions_mrmvar.R")

##load data
original_data=read.csv("molecular_list.csv")
group_data=read.csv("sample_list.csv")
names(group_data)[1]="id"
row.names(original_data)=original_data$formula

prokaryote=read.table("resample_16s.txt",sep="\t",head=1,row.names = 1)
output=original_data[colnames(prokaryote)]
input_ori=original_data[group_data$id[group_data$type=="Input"]]

input=output
for (ii in 1:42) {
  site=group_data$group[group_data$id==colnames(output)[ii]]
  outid=group_data$id[group_data$group==site&group_data$type=="Input"]
  input[ii]=input_ori[outid]
}

##match
prokaryote2=prokaryote
prokaryote3=prokaryote2[rowSums(prokaryote2)>0,]
prokaryote4=prokaryote3
for (i in 1:ncol(prokaryote4)){
  prokaryote4[i]=prokaryote4[i]/sum(prokaryote4[i])
}

input2=input[colnames(prokaryote4)]
input3=input2[rowSums(input2)>0,]

output2=output[colnames(prokaryote4)]
output3=output2[rowSums(output2)>0,]


##step MRM
###function
step_dis=function (data1,from1=100,to1=800,step1=1,halfweidth=20,method="bray",binary=F){
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

###set parameter
data1=output3
data2=original_data[c(1,7)]

from1=100
to1=800
step1=5
halfweidth=20
method="bray"
binary=F

prokary_d=vegdist(t(prokaryote4),method="bray",binary = F)
prokary_d_s=stdize(prokary_d)

input_d=vegdist(t(input3),method="bray",binary = F)
input_d_s=stdize(input_d)

x1=input_d_s
x2=prokary_d_s

data2=data2[row.names(data1),]

loopmain=seq((from1+halfweidth),(to1-halfweidth),step1)
loopsave=data.frame(loopmain,x1=0,x1x2=0,x2=0,unexplained=0)#save result
len=length(loopmain)

###calculate
for (ii in 1:len) {
  para1=loopmain[ii]-halfweidth
  para2=loopmain[ii]+halfweidth
  dothework=data1[data2[data2[2]>=para1&data2[2]<=para2,1],]#select data
  
  dothework_d=vegdist(t(dothework),method=method,binary = binary)
  dothework_d_s=stdize(dothework_d)
  
  mrmvar=varpart3(dothework_d_s,x1,x2)
  loopsave[ii,2:5]=mrmvar$part$indfract[,3]
}

outresult=loopsave

outresult[outresult<0]=0

outresult$unexplained=1-outresult$x1-outresult$x2-outresult$x1x2
