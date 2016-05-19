##recalculate ic50 vals

source("../../bin/ncatsSingleAgentScreens.R")
##goal is to collect data and format it into format that can be called by the add_gr_column.py script
require(dplyr)
require(tidyr)
createFile<-function(tab){
  tlist<-apply(tab,1,function(x){
    concs=x[26:36]
    counts=x[15:25]
   # zdata<-data.frame(cell_line=x[3],agent='-',perturbation=0,replicate=1,time=0,concentration=0,cell_count=x[8])
    cdata<-data.frame(cell_line=rep(x[3],length(concs)),
              agent=rep(x[37],length(concs)),
              perturbation=rep(0,length(concs)),
              replicate=rep(1,length(concs)),
              time=rep(48,length(concs)),
              concentration=concs,
              cell_count=counts,cell_count__ctrl=x[8],cell_count__time0=x[8])
    #return(rbind(zdata,cdata))
  })
  return(do.call('rbind',tlist))           
  
}

##take the allfiles object and combine into a single file
#slfiles<-lapply(allfiles,function(x) x[,c('Cell.line','name',)])
full.df<-do.call('rbind',lapply(allfiles,createFile))
write.table(full.df,'ncatsDataForgr50.tsv',sep='\t',quote=F,row.names=F)