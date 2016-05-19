source('../../bin/RNASeqData.R')

##all.ecounts<-cbind(ecounts[which(!is.na(t.idx)),],ecounts.mat[t.idx[which(!is.na(t.idx))],])
all.cors<-read.table(synGet('syn5595157')@filePath)

doCorPlots<-function(tcor,prefix=''){
  require(pheatmap)
    #first compute correlation
#  tcor=cor(df)
  pheatmap(tcor,main=paste(prefix,'Pearson Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'pearsonCors.pdf',sep='_'))
  
  #second: spearman
  tcor=cor(df,method='spearman')
  pheatmap(tcor,main=paste(prefix,'Spearman Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'spearmanCors.pdf',sep='_'))
  
  #third: quanile normalization
  
  tcor=cor(df,method='spearman')
  pheatmap(tcor,main=paste(prefix,'Spearman Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'spearmanCors.pdf',sep='_'))
  
  require(preprocessCore)
  qnormed=normalize.quantiles(as.matrix(df))
  colnames(qnormed)<-colnames(df)
  tcor=cor(qnormed)
  pheatmap(tcor,main=paste(prefix,'Q-normed Pearson Correlation'),cellheight=10,cellwidth=10,filename=paste(prefix,'qnormPearsonCor.pdf',sep='_'))
}



##add method to do some statistics on the correlations
doCorStats<-function(tcor,prefix='',minCor=0.25){
  #first compute correlation
 
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  
    doCorPlots(mmat,paste(prefix,'minPearsonCor',minCor,sep=''))
  
  #second: spearman
  tcor=cor(df,method='spearman')
  mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))
  sel.vals<-tcor[mvals,]
  
  write.table(sel.vals,paste(prefix,'spearmanCors.tsv',sep=''),sep='\t')
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minSpearmanCor',minCor,sep=''))
  
    #third: quanile normalization
  
 
  require(preprocessCore)
  qnormed=normalize.quantiles(as.matrix(df))
  colnames(qnormed)<-colnames(df)
  tcor=cor(qnormed)
  mvals=union(grep("ip",rownames(tcor)),grep('NF',rownames(tcor)))

  sel.vals<-tcor[mvals,]
  write.table(sel.vals,paste(prefix,'qnormedPearson.tsv',sep=''),sep='\t')
  mmat<-df[,colnames(sel.vals)[unique(which(sel.vals>minCor,arr.ind=T)[,2])]]
  doCorPlots(mmat,paste(prefix,'minQnormedPearsonCor',minCor,sep=''))              
  
}

##not looking great
doCorPlots(all.tpms,'tpm')
#doCorStats(all.ecounts,'estCounts')
