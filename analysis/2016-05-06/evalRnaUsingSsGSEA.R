##apply gene expression profiles to GSEA/MSIGDB using GSVA
library(GSVA)

source('../../bin/RNASeqData.R')


library(data.table)
library(GSEABase)
library(GSVAdata)
library(pheatmap)
#cut and paste from vignette
#cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,
#                            +                            min.sz=10, max.sz=500, verbose=TRUE)$es.obs,
#        +                            dir=cacheDir, prefix=cachePrefix)

data(c2BroadSets)

#> gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)$es.obs
rna.data<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)

#get entrez ids for RNA count matrix
eids<-as.data.frame(fread('../../data/HugoGIDsToEntrez_DAVID.txt',sep='\t',header=T))

eid.match=eids[match(rownames(rna.data),eids[,1]),2]
eid.rows=which(!is.na(eid.match))
ent.rna<-rna.data[eid.rows,]
rownames(ent.rna)=eid.match[eid.rows]

load('../../data/geneSets.RData')

es<-gsva(ent.rna,c2BroadSets,rnaseq=T,no.bootstraps=100,method='ssgsea')
eh<-gsva(rna.data,geneSet.hallmarks,rnaseq=T,no.bootstraps=100,method='ssgsea')
ec<-gsva(rna.data,geneSet.c6,rnaseq=T,no.bootstraps=100,method='ssgsea')
##let's write out the pathways to a tab-delimited file for future analysis


##pathway analysis directory
synid='syn5688828'
path.vals=es
#path.bootstrap=es$bootstrap$p.vals.sign
write.table(path.vals,'ssGSEAcellLineBroadPathEnrich.tab',sep='\t')

cpath.vals=eh
#path.bootstrap=es$bootstrap$p.vals.sign
write.table(cpath.vals,'ssGSEAcellLineHallmakrsPathEnrich.tab',sep='\t')

write.table(ec,'ssGSEAcellLineC6PathEnrich.tab',sep='\t')

#write.table(path.bootstrap,'cellLineBroadPath1000BootstrapPvals.tab',sep='\t')



#let's cluster, see how that goes? 
res=path.vals
plot(hclust(as.dist(1-cor(res)),method='ward.D2'))
colnames(res)[1]="ipNF05.5 (mixed clone)" 
colnames(res)[11]="ipNF05.5 (single clone)"

genotype=samp.mappings[match(colnames(res),samp.mappings[,1]),4]
names(genotype)<-colnames(res)

vars<-apply(res,1,var,na.rm=T)
mostvar<-res[order(vars,decreasing=T)[1:100],]
#names(pats)<-gsub('CT0*','Patient',rna.annot$synapseId)
pheatmap(mostvar,cellheight=10,cellwidth=10, annotation_col=data.frame(Genotype=genotype),
         filename='mostVariablePathways_ssGSEA.png')

c6res=ec
colnames(c6res)[1]="ipNF05.5 (mixed clone)" 
colnames(c6res)[11]="ipNF05.5 (single clone)"

vars<-apply(c6res,1,var,na.rm=T)
mostvar<-c6res[order(vars,decreasing=T)[1:50],]
#names(pats)<-gsub('CT0*','Patient',rna.annot$synapseId)
pheatmap(mostvar,cellheight=10,cellwidth=10, annotation_col=data.frame(Genotype=genotype),
         filename='mostVariableC6pathways_ssGSEA.png')

cres=cpath.vals
plot(hclust(as.dist(1-cor(cres)),method='ward.D2'))
colnames(res)[1]="ipNF05.5 (mixed clone)" 
colnames(res)[11]="ipNF05.5 (single clone)"

genotype=samp.mappings[match(colnames(cres),samp.mappings[,1]),4]
names(genotype)<-colnames(cres)

vars<-apply(cres,1,var,na.rm=T)
mostvar<-cres[order(vars,decreasing=T)[1:50],]
#names(pats)<-gsub('CT0*','Patient',rna.annot$synapseId)
pheatmap(mostvar,cellheight=10,cellwidth=10, annotation_col=data.frame(Genotype=genotype),
         filename='mostVariableHallmarkPathways_ssGSEA.png')



require(limma)

design= model.matrix(~factor(genotype))
colnames(design)=c('','+-','++')
fit <- lmFit(res, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="++")
sig=res[rownames(allGenes)[which(allGenes$P.Value<0.05)],]
pheatmap(sig,cellheight=10,cellwidth=10,annotation_col=data.frame(Genotype=genotype),
         filename='sigDiff_HomoZ_Pathways_ssGSEA.png')

fit <- lmFit(cres, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="++")
sig=cres[rownames(allGenes)[which(allGenes$P.Value<0.1)],]
pheatmap(sig,cellheight=10,cellwidth=10,annotation_col=data.frame(Genotype=genotype),
         filename='sigDiff_HomoZ_Hallmarks_ssGSEA.png')

fit <- lmFit(c6res, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="++")
sig=c6res[rownames(allGenes)[which(allGenes$P.Value<0.1)],]
pheatmap(sig,cellheight=10,cellwidth=10,annotation_col=data.frame(Genotype=genotype),
         filename='sigDiff_HomoZ_c6_ssGSEA.png')

flist=c('ssGSEAcellLineBroadPathEnrich.tab','ssGSEAcellLineHallmakrsPathEnrich.tab','sigDiff_HomoZ_Pathways_ssGSEA.png','sigDiff_HomoZ_Hallmarks_ssGSEA.png','mostVariablePathways_ssGSEA.png','mostVariableHallmarkPathways_ssGSEA.png')
flist=c('ssGSEAcellLineC6PathEnrich.tab','sigDiff_HomoZ_c6_ssGSEA.png','mostVariableC6pathways_ssGSEA.png')
for(file in flist){
  sf=File(file,parentId=synid)
  synStore(sf,activityName='ssGSEA enrichment analysis',
           used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-05-06/evalRnaUsingSsGSEA.R',wasExecuted=TRUE)))
}
