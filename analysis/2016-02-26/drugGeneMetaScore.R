##get PC values for each drug/cell combo
source("../../bin/ncatsSingleAgentScreens.R")

vars<-c("FAUC",'MAXR','LAC50','HILL','ZERO')
all.vals<-lapply(vars,getValueForAllCells)
names(all.vals)<-vars
drugs=rownames(all.vals$FAUC)
cells=colnames(all.vals$FAUC)

##now compute PCs and get variance

##how well do first pcs correlate between cell by drug vs. drug by cell?
cell.by.drug=sapply(cells,function(cell){
        drug.mat=sapply(drugs,function(drug){
           unlist(lapply(all.vals,function(x) x[[drug,cell]]))
          
       })
        drug.mat[which(is.na(drug.mat))]<-0.0
        pc=prcomp(t(drug.mat))
        pc$x[,1]        })

drug.by.cell<-sapply(drugs,function(drug){
  drug.mat=sapply(cells,function(cell){
    unlist(lapply(all.vals,function(x) x[[drug,cell]]))
    
  })
  drug.mat[which(is.na(drug.mat))]<-0.0
  pc=prcomp(t(drug.mat))
  pc$x[,1]        })

drug.cors=cor(drug.by.cell,t(cell.by.drug))
cell.cors=cor(t(drug.by.cell),cell.by.drug)

pdf("drugByCellEvaluation.pdf")
par(mfrow=c(1,2))
hist(drug.cors,main='Drug PC1 correlations')
hist(cell.cors,main='Cell PC1 correlations')
dev.off()

#we have to choose one dimension!!!  

##get RNA-Seq
source("../../bin/RNASeqData.R")
rna.dat<-rnaGencodeKallistoMatrix(useCellNames=TRUE,byGene=TRUE)
comm.cells<-intersect(cells,colnames(rna.dat))

r.dcor=cor(t(rna.dat[,comm.cells]),drug.by.cell[comm.cells,])
r.ccor=cor(t(rna.dat[,comm.cells]),t(cell.by.drug[,comm.cells]))

#now get targets
targs<-ncatsDrugTargets()

targ.df=targs
targ.df$RnaDrugCor=apply(targs,1,function(x) 
  if(x[[1]]%in%colnames(r.dcor)&&x[[2]]%in%rownames(r.dcor)) 
    return(r.dcor[x[[2]],x[[1]]]))

targ.df$RnaCellCor=apply(targs,1,function(x) 
  if(x[[1]]%in%colnames(r.ccor)&&x[[2]]%in%rownames(r.ccor)) 
    return(r.ccor[x[[2]],x[[1]]]))

