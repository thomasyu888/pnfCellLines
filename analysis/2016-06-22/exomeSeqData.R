source('../../bin/exomeSeqProcessing.R')

all.dat<-getAnnovarFiles(reStore=TRUE)

nf1.muts<-subset(all.dat,KnownGeneGeneName=="NF1")

library(ggplot2)

ggplot(nf1.muts)+geom_point(aes(x=StartPosition,y=CellLine,color=KnownGeneGeneLocation,shape=KnownGeneExonFunction))
ggsave('nf1Mutants.png')

##now try to collect more detailed mutants.
ex.muts<-subset(all.dat,KnownGeneGeneLocation=='exonic')
ns.muts<-subset(ex.muts,KnownGeneExonFunction!='synonymous SNV')
het.muts<-subset(ns.muts,Genotype=='het')
hom.muts<-subset(ns.muts,Genotype=='hom')

require(dplyr)
require(reshape2)
het.counts<-het.muts%>%group_by(CellLine)%>%summarize(numGenes=n_distinct(KnownGeneGeneName))

hom.counts<-hom.muts%>%group_by(CellLine)%>%summarize(numGenes=n_distinct(KnownGeneGeneName))

##now create table
hom.mat<-acast(hom.muts,KnownGeneGeneName~CellLine,value.var="EnsgeneExonFunction")
het.mat<-acast(het.muts,KnownGeneGeneName~CellLine,value.var="EnsgeneExonFunction")

hom.shared=which(apply(hom.mat,1,function(x) all(x>0)))
het.shared=which(apply(het.mat,1,function(x) all(x>0)))

png('exomeSeqClustering.png',height=800)
par(mfrow=c(2,1))
plot(hclust(dist(t(het.mat)),method='ward'),main='Heterozygous non-synonymous variants')
plot(hclust(dist(t(hom.mat)),method='ward'),main='Homozygous non-synonymous variants')
dev.off()

synStore(File('exomeSeqClustering.png',parentId='syn6174634'),
         executed=list(list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2016-06-22/exomeSeqData.R')),
                used=list(list(entity='syn6086887')))