##mor drug analysis


source('../../bin/crossDataComps.R')

##first cluster all TAUC values across cell lines, plot in heatmap
tauc=getValueForAllCells("TAUC")



fauc=getValueForAllCells('FAUC')
nzfauc=fauc
nzfauc[which(is.na(nzfauc),arr.ind=T)]<-1000.0
#pheatmap(nzfauc,cellheight=10,cellwidth=10,file='FAUC_for_all_cells.png')
#pheatmap(nzfauc,cellheight=10,cellwidth=10,file='FAUC_for_all_cells.pdf')

zsfauc<-apply(fauc,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
nz.zs.fauc=zsfauc
nz.zs.fauc[which(is.na(nz.zs.fauc),arr.ind=T)]<-0.0
#pheatmap(nz.zs.fauc,cellheight=10,cellwidth=10,file='Zscore_FAUC_for_all_cells.pdf')

##there is a pretty blue cluster there, can we do any enrichment? 
drug.dists<-dist(nz.zs.fauc)
hc=hclust(drug.dists)


require(dplyr)
require(reshape2)

drug.clusters<-cutree(hc,h=4)
df=data.frame(Drug=names(drug.clusters),Cluster=as.numeric(drug.clusters),Target=drugTargs[names(drug.clusters)])
ndf<-subset(df,!is.na(Target))


##create data frame calculating all cluster statistics for each drug target
drug.targs= ndf %>% group_by(Target,Cluster) %>% summarise(TimesTargetInCluster=n())
tot.targs<-ndf %>% group_by(Target) %>% summarize(Total=n())
clust.size<-ndf %>% group_by(Cluster) %>% summarize(Total=n())

drug.targs$NumDrugsWithTarget=tot.targs$Total[match(drug.targs$Target,tot.targs$Target)]
drug.targs$ClusterSize=clust.size$Total[match(drug.targs$Cluster,clust.size$Cluster)]
targ.size=drug.targs %>% group_by(Cluster)%>% summarize(NumTargets=n())
drug.targs$UniqueTargets=targ.size$NumTargets[match(drug.targs$Cluster,targ.size$Cluster)]

drug.targs$Drug=ndf$Drug[match(drug.targs$Cluster,ndf$Cluster)]
pvals=apply(drug.targs,1,function(x){
  over=as.numeric(x[['TimesTargetInCluster']])
  tt=as.numeric(x[['NumDrugsWithTarget']])
  cs=as.numeric(x[['ClusterSize']])
 mat=matrix(c(over,tt-over,cs-over,nrow(ndf)-tt-cs+over),nrow=2)
 return(fisher.test(mat,alt='g')$p.value)
})

drug.targs$Pvalue=pvals
drug.targs$FDR=p.adjust(pvals,method='fdr')

pmat<-acast(drug.targs,Target~Cluster,value.var="Pvalue",fill=1.0)
qmat<-acast(drug.targs,Target~Cluster,value.var="FDR",fill=1.0)

#colnames(pvals)=as.character(clust.size$Cluster)
#rownames(pvals)<-as.character(tot.targs$Target)
#drug.targs$Pvalues=apply(drug.targs,1,function(x) pvals[match(x[[1]],rownames(pvals)),match(x[[2]],colnames(pvals))])



pheatmap(log10(pmat),file='log10PvaluesOfClusters.png',cellheight=10,cellwidth=10)
pheatmap(log10(qmat),file='log10FDROfClusters.png',cellheight=10,cellwidth=10)

write.table(drug.targs,'drugTargetClustersAndEnrichment.txt',col.names=T,row.names=F,sep='\t')

