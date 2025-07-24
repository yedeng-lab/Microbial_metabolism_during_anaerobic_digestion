#Counts of shared molecules in digesters
#@ yang
rm(list=ls())

##function
group_mean=function(data,group){
  name=as.data.frame(table(group[1]))
  out=data[,1:nrow(name)];colnames(out)=name[,1]
  for(i in 1:nrow(name)){
    list=!is.na(match(group[,1],name[i,1]))
    out[,i]=rowMeans(data[,match(row.names(group)[list],colnames(data))])
  }
  return(out)
}

##load data
original_data=read.csv("molecular_list.csv")
group_data=read.csv("sample_list.csv")
names(group_data)[1]="id"
row.names(group_data)=group_data$id
row.names(original_data)=original_data$formul

## data of reactors
data=original_data[c(group_data$id[group_data$type=="Output"])]
data=data[rowSums(data)>0,]
data[data!=0]=1
group2=group_data[group_data$type=="Output",2:3]
group2=group2[-2]

mean_count=group_mean(data,group2)

result1=as.data.frame(matrix(NA,ncol=3,nrow=3))
colnames(result1)=c("Group","Type","Number")
result1$Group="All"
result1$Type=c("Occurring","Ubiquitous","Persistent")
result1$Number[1]=nrow(data)
mean_count2=mean_count
mean_count2[mean_count!=0]=1
result1$Number[2]=nrow(mean_count2[rowSums(mean_count2)>=7,])
result1$Number[3]=nrow(mean_count[rowSums(mean_count)>=6.9,])

for (ii in 1:7) {
  result2=as.data.frame(matrix(NA,ncol=3,nrow=2))
  colnames(result2)=c("Group","Type","Number")
  result2$Type=c("Occurring","Persistent")
  result2$Group=names(mean_count)[ii]
  result2$Number[1]=nrow(mean_count[mean_count[ii]>0,])
  result2$Number[2]=nrow(mean_count[mean_count[ii]>0.9,])
  
  result1=rbind(result1,result2) 
}
