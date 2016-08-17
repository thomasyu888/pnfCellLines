##reformat gene matrix to compile by gene and sample name
source('../../bin/RNASeqData.R')
mat<-rnaGencodeKallistoMatrix(byGene=TRUE,useCellNames=T)