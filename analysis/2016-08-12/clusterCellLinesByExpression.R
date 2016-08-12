#' Cluster cell lines by tissue type by expression
#'
source("../../../RASPathwaySig/bin/cBioPortalData.R")

##get zscored data
ccle.mat<-getCcleExpressionData('',getZscores=TRUE)

#map by cell line name
all.tiss<-unique(sapply(colnames(ccle.mat),function(x) paste(unlist(strsplit(x,split='_'))[-1],collapse='_')))
ccle.tiss.averages<-sapply(all.tiss,function(x) {
  cols=grep(x,colnames(ccle.mat))
  if(length(cols)>1)
    return(rowMeans(ccle.mat[,cols],na.rm=T))
  else
    return(ccle.mat[,cols])})

##get pnf cell line data
require(synapseClient)
synapseLogin()

#get the TPM values
tpm.vals<-read.table(synGet('syn5580347')@filePath,sep='\t',header=T)
#break down by gene

genes<-unique(sapply(rownames(tpm.vals),function(x) unlist(strsplit(x,split='.',fixed=T))[1]))

per.gene.vals<-t(sapply(genes,function(x){
  vals<-grep(paste(x,'.',sep=''),rownames(tpm.vals),fixed=T)
  if(length(vals)>1)
    return(colSums(tpm.vals[vals,],na.rm=T))
  else
    return(tpm.vals[vals,])
  
}))

zscore<-function(x){
  x<-unlist(x)
  (x-mean(x,na.rm=T))/sd(x)
}
#z-score by gene
z.score.gene<-apply(per.gene.vals,2,zscore)

##get column names
cell.line.names<-synTableQuery('SELECT "Sample Name","RNA-Seq Data (Gencode)" FROM syn5014742 where "RNA-Seq Data (Gencode)" is not NULL')@values
idx<-match(colnames(z.score.gene),cell.line.names$`RNA-Seq Data (Gencode)`)

colnames(z.score.gene)<-cell.line.names$`Sample Name`[idx]

##get other nf data
