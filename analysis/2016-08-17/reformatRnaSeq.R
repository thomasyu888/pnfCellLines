##reformat gene matrix to compile by gene and sample name
source('../../bin/RNASeqData.R')
mat<-rnaGencodeKallistoMatrix(byGene=TRUE,useCellNames=T)

write.table(mat,file='pnfGencodeKallisto_tpmByGene.tsv',sep='\t',row.names=T,col.names=T)
synStore(File('pnfGencodeKallisto_tpmByGene.tsv',parentId='syn5579785'),executed='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2016-08-17/reformatRnaSeq.R',used='syn5580347')

##now create an expression set. 
library(Biobase)
expr<-ExpressionSet(mat)

annotes=synTableQuery('SELECT "Sample Name","Sample Genotype","RNA-Seq Data (Gencode)" FROM syn5014742 where "RNA-Seq Data (Gencode)" is not null')@values

df<-data.frame(annotes[,2:3])
rownames(df)<-annotes[,1]
colnames(df)<-c('Genotype','Synapse ID')
pData(expr)<-df
write.table(df,row.names=T,col.names=T,sep='\t',file='pnfPhenoData.tsv')

synStore(File('pnfPhenoData.tsv',parentId='syn5579785'),executed='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2016-08-17/reformatRnaSeq.R',used='syn5580347')
