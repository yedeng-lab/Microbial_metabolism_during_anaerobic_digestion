# random forest
# to identify molecules in dIfferent group
# input VS. reactor
# @ yang

rm(list=ls())

##load data
original_data=read.csv("molecular_list.csv")
group_data=read.csv("sample_list.csv")
names(group_data)[1]="id"
purple_data=original_data[group_data$id]
use_group=group_data[3]
rownames(use_group)=group_data$id

tax_data=original_data[c("formula","classification")];row.names(tax_data)=original_data$formula

## prepation for randomForest
library(microeco)
library(randomForest)
dataset=microtable$new(sample_table = use_group,
                       otu_table = purple_data,
                       tax_table = tax_data)

##run randomtree
set.seed(8080)
rf=trans_diff$new(dataset = dataset,
                  method = "rf",
                  group="type",
                  taxa_level = "formula")

outfile=rf$res_diff