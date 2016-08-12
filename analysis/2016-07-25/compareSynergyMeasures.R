##compare synergy across measures
source("../../bin/getSynergyAndStats.R")
source("../../bin/ncatsCombinationScreens.R")
library(ggplot2)
for(combo in c('6x6','10x10')){
  cells<-getCellForCombo(combo)
  measure='CTG'
  for(cell in cells[,1]){
    ncats<-getNcatsSynergyScores(cell,combo,measure)
    recalc<-reCalculateSynergyScores(cell,combo,measure)

    ##first combine everything into a single data frame
    full.tab<-ncats%>%select(match(c("RowName","ColName","Method","Value"),colnames(ncats)))%>%bind_rows(recalc)
    #get the drug combos as a unque column
    new.tab<-mutate(full.tab,Combo=paste(RowName,ColName,sep='_'))     
    method.cors<-cor(acast(new.tab,Combo~Method,value.var="Value",fun.aggregate=mean),use='pairwise.complete.obs')
    pheatmap(method.cors,cellwidth = 10,cellheight=10,main=paste('Method correlations in',cell),filename=paste('MethodCorrelationsIn',gsub('\\)|\\(| ','',cell),'.png',sep=''))
    
    }
}