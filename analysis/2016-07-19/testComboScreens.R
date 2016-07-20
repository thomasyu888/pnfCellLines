source("../../bin/ncatsCombinationScreens.R")

testdata<-getResponseData('6x6','CTG')
metadata<-getMetadataForCombo('6x6','CTG')
syndata<-getSynergyData('6x6','CTG')

##got some data, for the pacakge we need the responses
library(synergyfinder)

tmat<-testdata$HFF$`1`

##looks like we need to order the rows and columns
otmat<-tmat[order(as.numeric(rownames(tmat))),order(as.numeric(colnames(tmat)))]

kmat<-100-otmat

