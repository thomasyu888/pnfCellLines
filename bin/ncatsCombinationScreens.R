##begin to evaluate combination screens....

##
## Drug sensitivity data files for NCATS single agent screens
##
##
library(synapseClient)
library(data.table)
library(pheatmap)
synapseLogin()
library(ggplot2)
require(dplyr)
require(tidyr)
require(reshape2)
screendirs=list(`10x10`='syn5611797',`6x6`='syn5611796')


#'Get list of all cells available for combo and measure
#'@param comboScreen
#'@return data frame from synapse table query
getCellForCombo<-function(comboScreen=c('10x10','6x6'))
{
  fileparent=screendirs[[comboScreen]]  
  cells=synapseQuery(paste("select name,id from entity where parentId=='",fileparent,"'",sep=''))
  print(paste('Retrieved',nrow(cells),'cell types for',comboScreen,'screen'))
  return(cells)
}

#'Get associated metadata for cell, combo and measure
#'@param cell name of cell to be evaluated
#'@param comboScreen either 6x6 or 10x10
#'@param measure either CTG or CCG (for 10x10 only)
#'@return data frame of stats for particular cell 
getMetadataForCombo<-function(cell='HFF',comboScreen=c('10x10','6x6'),measure=c('CTG','CCG')){
    if(measure=='CCG'&&comboScreen=='6x6'){
      print('CCG data not available for 6x6, evaluating with CTG')
      measure='CTG'
    }
      
#    cell.dat<-lapply(cells[,1],function(x){
    cells<-getCellForCombo(comboScreen)
      synid=cells[match(cell,cells[,1]),2]
     # print(synid)
      allf<-synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
      allf<-allf[grep(measure,allf[,1]),]
      meta<-data.frame(fread(synGet(allf[grep('metadata',allf[,1]),2])@filePath))
      calc<-data.frame(fread(synGet(allf[grep('calc',allf[,1]),2])@filePath))
      ##not sure which values to use?  
      combined<-cbind(meta,calc)
      ##return rowtrug, row target, col drug, col target, and some synergy score....
      return(combined)
 
}


#' Get response matrix for combination screen
#' 
#' @param cell - name of cell to be evaluated
#' @param comboScreen - either 6x6 or 10x10
#' @param measure - either CTG or CCG (for 10x10only)
#' @return matrix representing the fraction of the cells remaining
getResponseData<-function(cell='HFF',comboScreen=c('10x10','6x6'),measure=c('CTG','CCG')){
  if(measure=='CCG'&&comboScreen=='6x6'){
    print('CCG data not available for 6x6, evaluating with CTG')
    measure='CTG'
  }

  cells<-getCellForCombo(comboScreen)
  
 # cell.dat<-lapply(cells[,1],function(x){
    synid=cells[match(cell,cells[,1]),2]
    allf<-synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
    allf<-allf[grep(measure,allf[,1]),]
    
    meta<-data.frame(fread(synGet(allf[grep('metadata',allf[,1]),2])@filePath))
    
    resp<-data.frame(fread(synGet(allf[grep('resp',allf[,1]),2])@filePath))
    ##not sure which values to use?  
    all.mats<-lapply(unique(resp$BlockId),function(x) acast(subset(resp,BlockId==x),Row~Col,value.var="Value",fun.aggregate=mean))
    names(all.mats)<-unique(resp$BlockId)
    
    named.mats<-lapply(names(all.mats),function(x){
      mval<-all.mats[[x]]
      metas<-subset(meta,BlockId==x)
      cols<-unlist(strsplit(metas$ColConcs,split=','))
      rows<-unlist(strsplit(metas$RowConcs,split=','))
      rownames(mval)<-rows
      colnames(mval)<-cols
      return(mval)
    })
    names(named.mats)<-names(all.mats)
    
    ##return rowtrug, row target, col drug, col target, and some synergy score....
    return(named.mats)
  
}


#'Format synergy data to be used by the synergyfinder R package
#'There is a bug in the synergy calculators so try to use https://github.com/sgosline/synergyfinder
#'@param cell cell of interest
#'@param comboScreen either 6x6 or 10x10
#'@param measure either CTG or CCG (for 10x10 only)
#'@return list of data to be analyzed by synergyfinder
#'
getSynFinderData<-function(cell='HFF',comboScreen=c('10x10','6x6'),measure=c('CTG','CCG')){
  if(measure=='CCG'&&comboScreen=='6x6'){
    print('CCG data not available for 6x6, evaluating with CTG')
    measure='CTG'
  }
  

#  cell.dat<-lapply(cells[,1],function(x){
  
  synid=cells[match(cell,cells[,1]),2]
  allf<-synapseQuery(paste("select name,id from entity where parentId=='",synid,"'",sep=''))
  allf<-allf[grep(measure,allf[,1]),]
    
  meta<-data.frame(fread(synGet(allf[grep('metadata',allf[,1]),2])@filePath))
    
  resp<-data.frame(fread(synGet(allf[grep('resp',allf[,1]),2])@filePath))

    ##first separate out row and column concentrations
    if(comboScreen=='6x6'){
      rowCon=paste('row',c(1,2,3,4,5,6),sep='_')
      colCon=paste('col',c(1,2,3,4,5,6),sep='_')
    }else{
      rowCon=paste('row',c(1,2,3,4,5,6,7,8,9,10),sep='_')
      colCon=paste('col',c(1,2,3,4,5,6,7,8,9,10),sep='_')
    }
    fullDf<-tidyr::separate(meta,RowConcs,into=rowCon,sep=',')%>%tidyr::separate(ColConcs,into=colCon,sep=',')
    
    if(comboScreen=='6x6'){
      rowdf<-fullDf%>%select(match(c("BlockId",rowCon),names(fullDf)))%>%tidyr::gather('Row','RowConc',2:7)%>%arrange(BlockId)
      coldf<-fullDf%>%select(match(c("BlockId",colCon),names(fullDf)))%>%tidyr::gather('Col','ColConc',2:7)%>%arrange(BlockId)
    } else{
      rowdf<-fullDf%>%select(match(c("BlockId",rowCon),names(fullDf)))%>%tidyr::gather('Row','RowConc',2:11)%>%arrange(BlockId)
      coldf<-fullDf%>%select(match(c("BlockId",colCon),names(fullDf)))%>%tidyr::gather('Col','ColConc',2:11)%>%arrange(BlockId)
      
    }   
    ##now rename rows and columns
    rowdf$Row<-sapply(rowdf$Row,function(x) gsub('row_','',x))
    coldf$Col<-sapply(coldf$Col,function(x) gsub('col_','',x))
    
    metdf<-select(fullDf,match(c('BlockId','RowName','ColName','RowConcUnit','ColConcUnit'),colnames(fullDf)))
    
    ##join rows and columns, sort both the metadata and responses in the same order
    jm<-full_join(metdf,rowdf,by='BlockId')%>%full_join(coldf,by='BlockId')%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
    resp<-resp%>%arrange(Col)%>%arrange(Row)%>%arrange(BlockId)
    
    final.df<-bind_cols(jm,select(resp,c(Value,Replicate)))%>%select(c(BlockId,Row,Col,Value,Replicate,RowName,ColName,RowConc,ColConc,RowConcUnit))
    colnames(final.df)<-c('BlockID','Row','Col','Response','Replicate','DrugRow','DrugCol','ConcRow','ConcCol','ConcUnit')
    final.df$ConcRow<-as.numeric(final.df$ConcRow)
    final.df$ConcCol<-as.numeric(final.df$ConcCol)
    ##return rowtrug, row target, col drug, col target, and some synergy score....
    return(final.df)

}



#'Plot synergy metdata values across cells - experiment with the best way to visualize the data in its 
#'@param file.list list of results from getMetadataForCombo function, nam
#'@param prefix to add to file
#'@param value Which synergy value to collect from metadata
#'@return data frame of all results
plotValsAcrossCells<-function(file.list,prefix='',value='Beta'){
  fres<-do.call('rbind',lapply(names(file.list),function(y){
    x<-file.list[[y]]
    data.frame(Row=x$RowName,Col=x$ColName,
            RowTarg=x$RowTarget,ColTarg=x$ColTarget,
            Value=x[,value],Cell=y)
  }))
  
  
  fname=paste(prefix,value,'Values',sep='')
    
  ##plot values by cell
  ggplot(fres)+geom_boxplot(aes(y=Value,x=Cell))
  ggsave(paste(fname,'ByCell.png',sep=''))
  
  ##plot values by row, column to see if there are outliers 
  ggplot(fres)+geom_boxplot(aes(y=Value,x=Row,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByRow.png',sep=''))
  
  ggplot(fres)+geom_boxplot(aes(y=Value,x=Col,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByColumn.png',sep=''))
  
  ##then row/column targets
  
  ggplot(fres)+geom_boxplot(aes(y=Value,x=RowTarg,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByRowTarget.png',sep=''))
  
  ggplot(fres)+geom_boxplot(aes(y=Value,x=ColTarg,fill=Cell))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(fname,'ByColTarget.png',sep=''))
  return(fres)
}


#'Performs a linear model to evaluate if there were speicfic combinations that were more synergistic
#'than others based on value assessed
#'@param value.df DataFrame returned by plotValsAcrossCells
#'@param prefix for file to save
#'@return list of drug coefficients and target coefficients
doLinearModel<-function(value.df,prefix){
  
  lres<-lm(Value~Row+Col+Cell+Row*Col,value.df)
  dcoeffs<-data.frame(summary(lres)$coefficients)
  #print(head(coeffs))
  dcoeffs$adjustedP<-p.adjust(as.numeric(dcoeffs[,4]))
  dcoeffs<-dcoeffs[order(dcoeffs[,4],decreasing=F),]
  write.table(dcoeffs,paste(prefix,'DrugLinearModelCoefficients.tsv',sep=''),sep='\t',row.names=T,col.names=T)
  
  lres<-lm(Value~RowTarg+ColTarg+Cell+RowTarg*ColTarg,value.df)
  tcoeffs<-data.frame(summary(lres)$coefficients)
  tcoeffs$adjustedP<-p.adjust(as.numeric(tcoeffs[,4]))
  tcoeffs<-tcoeffs[order(tcoeffs[,4],decreasing=F),]
  write.table(tcoeffs,paste(prefix,'DrugTargetLinearModelCoefficients.tsv',sep=''),sep='\t',row.names=T,col.names=T)
  
  return(list(drug=dcoeffs,target=tcoeffs))
  
}

#'Plot the results of the linear model
#'@param lmres results from doLinearModel
#'@param value.df Original input into same function
#'@param pval pvalue cutoff
#'
plotLMResults<-function(lmres,value.df,pval=0.05){
    combos<-lmres[grep(':',rownames(lmres)),]
    combos<-cbind(combos,t(sapply(rownames(combos),function(x) unlist(strsplit(gsub("Row|Col",'',x),split=':')))))
    colnames(combos)[6:7]<-c('Row','Col')    
    #pmat<-reshape2::dcast(combos,Row~Col,value.var='Pr...t..')
    apmat<-reshape2::dcast(combos,Row~Col,value.var='adjustedP')
    rownames(apmat)<-apmat$Row
    apmat<-apmat[,-1]
    apmat[which(is.na(apmat),arr.ind=T)]<-1.0
    
    mat<-reshape2::dcast(value.df,Row~Col,value.var='Value',fun.aggregate=mean)
    rownames(mat)<-mat$Row
    mat<-mat[,-1]
    mat[which(is.na(mat),arr.ind=T)]<-0.0
    rord<-order(apply(mat,1,sum))
    cord<-order(apply(mat,2,sum))
    
    par(mfrow=c(2,1))
    pheatmap(mat[rord,cord],cellheight=10,cellwidth=10,cluster_rows = F,cluster_cols = F)
}
    
    
    
    