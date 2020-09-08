rm(list=ls(all=TRUE))

library(e1071)
library(mclust)

###### SET YOUR FOLDER PATH #####
setwd="/ufrc/boucher/kingdgp/anthony-omsim/OMSIM_MODIFIED-omsim_modified/OMSIM_MODIFIED-omsim_modified/omsim/test/ecoli/myclustering/"
fileName="all_6mers.txt"

data=read.table(fileName,header=F,sep=" ",na="?")

data=data[,-dim(data)[2]]
names(data)=paste("kmer",1:dim(data)[2],sep="_")

###### ASINH(.) TRANSFORM #####
data=asinh(data)



###### FIT MODEL #####
bestn=1000
model=Mclust(data,G=bestn,modelNames="EII")
summ=summary(model, parameters=T)

###### SAVE MODEL PARAMETERS #####
capture.output(summ,file="model_summary.txt")

###### PRINT MODEL PREDICTIONS #####
pred=predict(model)
write.csv(pred$classification,file="model_predictions.csv")
#NOTE: IF YOU WANT PROBABILITY MEMBERSHIP TOO, USE write.csv(pred,file="model_predprob.csv")




