###a set of functions to get synergy using 


#'NCATS provides various metrics of synergy for each drug combo
#'@param cell
#'@param combo
#'@param measure
getNcatsSynergyScores<-function(cell,combo,measure){
  md<-getMetadataForCombo(cell,combo,measure)
  md<-md[,c('RowName','RowTarget','ColName','ColTarget','ExcessHSA','ExcessCRX','LS3x3','Beta','Gamma',"DBSumPos","DBSumNeg")]
  md<-gather(md,"Method","Value",5:11)
  return(md)
  
  }


#'We can also recalculate the synergy using the synergyfinder package.
#'Warning, this can take a while, so maybe pre-compute everything and store on synapse for the future.
#'@param cell
#'@param combo
#'@param measure
reCalculateSynergyScores<-function(cell,combo,measure){
  require(synergyfinder)
  dat<-getSynFinderData(cell,combo,measure)
  ##now we have to recalculate the three metrics.
  reshaped<-ReshapeData(dat)
  
  full.df<-do.call("rbind",lapply(c("ZIP", "HSA", "Bliss", "Loewe"),function(method){
    print(paste('Calculating synergy scores using',method,'algorithm...'))
    scores<-CalculateSynergy(reshaped,method=method)
    
    res<-cbind(scores$drug.pairs[,c('drug.row','drug.col')],Value=unlist(lapply(scores$scores,mean)),Method=rep(paste('mean',method,sep=''),length(scores$scores)))
    colnames(res)[1:2]<-c('RowName','ColName')
    return(res)
  }))
  return(full.df)
    
}

