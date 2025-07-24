#calculate transformation use moteam method
#@ yang 
rm(list=ls())

##load original data
original_data=read.csv("molecular_list.csv")
group_data=read.csv("sample_list.csv")
names(group_data)[1]="id"
purple_data=original_data[group_data$id]
row.names(purple_data)=original_data$formula

##load function
source("moteam_function.R")
source("get_dir_function.R")

##clean data
use_list=group_data[group_data$type=="Output",]
group_list=unique(group_data$group)

##main loop

for (ii in 1:7){
  ###get data by group
  group_name=group_list[ii]
  sample_name=use_list$id[use_list$group==group_name]
  thedata=purple_data[sample_name]
  print(names(thedata))
  ###get matrix
  min_occur=200
  p_adj_method="none"
  select_num=5
  threshold=0.99
  method="mic"
  used_p=0.05
  
  cal_moteam=moteam(thedata,n.majority=select_num,threshold=threshold,min_occur=min_occur,p_adj_method=p_adj_method,used_p=used_p,method = method)
  
  write.table(cal_moteam$MTE,paste0(group_name,"_MTE_",method,threshold,min_occur,p_adj_method,".txt"),sep="\t",quote = F,row.names = F)
  write.table(cal_moteam$network_matrix,paste0(group_name,"_network_",method,threshold,min_occur,p_adj_method,".txt"),sep="\t",quote = F,col.names =NA)
  write.table(cal_moteam$relation_table,paste0(group_name,"_reaction_",method,threshold,min_occur,p_adj_method,".txt"),sep="\t",quote = F,row.names = F)
  ###get direction
  link=cal_moteam$relation_table
  data1=thedata
  
  method="mic"
  n_expansion=1000000#The calculation of mic requires that the value is not too low, and that the multiplier does not affect the result.
  result=get_dir(link,data1,method,n_expansion)
  write.table(result,paste0(group_name,"_direction_",method,threshold,min_occur,p_adj_method,".txt"),sep="\t",quote = F,row.names = F)
}
