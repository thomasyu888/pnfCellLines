##do some basic drug sensitivity analysis

source("../../bin/drugSensData.R")

allcells=dfiles$entity.sampleName
##let's start with pretty plots of individual cells and response curve classes...
res<-lapply(allcells,plotOneCell,as.categ=TRUE)

res<-lapply(allcells,plotOneCell,use.disc=TRUE)

##now take each of the elements going into the calculations above and plot them across cell lines
#res<-lapply(valsOfInterest,plotMostVariableVals)

dv<-plotMostVariableVals('CRC','pdf')