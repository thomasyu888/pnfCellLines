##Get identifiers from exome seq data and map them to appropriate patient labels

library(synapseClient)
synapseLogin()

tab<-read.table(synGet('syn6092341')@filePath,sep=',',header=T)
