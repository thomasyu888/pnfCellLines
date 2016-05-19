##mor drug analysis


source('../../bin/crossDataComps.R')

##first cluster all TAUC values across cell lines, plot in heatmap
tauc=getValueForAllCells("TAUC")

pheatmap(tauc,cellheight=10,cellwidth=10,file='TAUC_for_all_cells.png')
pheatmap(tauc,cellheight=10,cellwidth=10,file='TAUC_for_all_cells.pdf')


fauc=getValueForAllCells('FAUC')
nzfauc=fauc
nzfauc[which(is.na(nzfauc),arr.ind=T)]<-1000.0
pheatmap(nzfauc,cellheight=10,cellwidth=10,file='FAUC_for_all_cells.png')
pheatmap(nzfauc,cellheight=10,cellwidth=10,file='FAUC_for_all_cells.pdf')

zsfauc<-apply(fauc,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
nz.zs.fauc=zsfauc
nz.zs.fauc[which(is.na(nz.zs.fauc),arr.ind=T)]<-0.0
pheatmap(nz.zs.fauc,cellheight=10,cellwidth=10,file='Zscore_FAUC_for_all_cells.pdf')

##there is a pretty blue cluster there, can we do any enrichment? 
drug.dists<-dist(nz.zs.fauc)
hc=hclust(drug.dists)

drug.clusters<-sapply(unique(cutree(hc,h=4)),function(x) names(which(cutree(hc,h=4)==x)))
hist(sapply(drug.clusters,length))

##do target enrichment


##now move onto RNA-Seq/correlation
##experiment with z-normalization of the pearson correlation values
drug.to.targ<-getValueForAllCells('target')[,1]
alltargs<-setdiff(unique(unlist(sapply(drug.to.targ,strsplit,', '))),NA)



