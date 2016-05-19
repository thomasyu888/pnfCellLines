###upload/process CCLE genetic data into meaningful format

library(synapseClient)
library(pheatmap)
synapseLogin()
source("../../bin/ctrpSingleAgentScreens.R")

storeMSSMprocessedMutationCalls<-function(){
  pid='syn5706496'
  fname='../2016-03-03/gene_attribute_matrix.txt'
  synStore(File(fname,parentId=pid))
}



#storeMSSMprocessedMutationCalls()

##now get nF1 mutations and compute differences between them 
ccle.calls=getMSSMprocessedMutationCalls()
source('../../bin/singleDrugAnalysis.R')

auc.mat<-getCtrpScreensAsMatrix()

overlap=intersect(colnames(auc.mat),colnames(ccle.calls))
c.calls=ccle.calls[,overlap]
gt=rep('+',length(overlap))
gt[which(c.calls['NF1',]==1)]<-'-'
names(gt)<-overlap

auc.mat[is.nan(auc.mat)]<-NA
auc.mat[is.na(auc.mat)]<-mean(auc.mat,na.rm=T)

tab<-aucDifferentialResponse(auc.mat[,overlap],gt)

sigs=rownames(tab)[which(tab$P.Value<0.005)]

#drug.targs<-ctrpDrugTargets()
#all.targs<-as.character(drug.targs$Target)
#names(all.targs)<-drug.targs$Drug

#write to file
write.table(tab,'aucValuesChangingAcrossNF1Status.tab')
pheatmap(auc.mat[sigs,overlap],#annotation_row=data.frame(Target=all.targs),
         clustering_distance_rows='correlation',
         annotation_col=data.frame(Genotype=gt),cellwidth=10,cellheight=10,file='aucValuesDistinctAcrossNF1Status_p005.png')



##now do that for the nplr-normalized data
auc.mat=ctrpDoseResponseCurve(recalculate=F,as.matrix=T)
overlap=intersect(colnames(auc.mat),colnames(ccle.calls))
gt=rep('+',length(overlap))
gt[which(ccle.calls['NF1',overlap]==1)]<-'-'
names(gt)<-overlap


drug.targs<-ctrpDrugTargets()
all.targs<-as.character(drug.targs$Target)
names(all.targs)<-drug.targs$Drug
auc.mat[is.nan(auc.mat)]<-NA
auc.mat[is.na(auc.mat)]<-mean(auc.mat,na.rm=T)
tab<-aucDifferentialResponse(auc.mat[,overlap],gt)

sigs=rownames(tab)[which(tab$P.Value<0.005)]

#write to file
write.table(tab,'recalculatedAucValuesChangingAcrossNF1Status.tab')
pheatmap(auc.mat[sigs,overlap],#annotation_row=data.frame(Target=all.targs[overlap]),
         annotation_col=data.frame(Genotype=gt[overlap]),#clustering_distance_rows = 'correlation',
         #clustering_distance_cols = 'correlation',
         cellwidth=10,cellheight=10,file='recalculatedAucValuesDistinctAcrossNF1Status_p005.png')

#lastly for the combined dataset. 