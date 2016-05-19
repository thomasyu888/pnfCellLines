##now we want to test the clustering parameters

source("../../bin/singleDrugAnalysis.R")
source("../../bin/ctrpSingleAgentScreens.R")
source("../../bin/ncatsSingleAgentScreens.R")



fs=synapseQuery("select name,id from entity where parentId=='syn5674273'")
fs=fs[grep('csv',fs[,1]),]
for(sid in fs[,2])
  synDelete(sid)

afiles=list.files('../2016-02-17')
cluster.dir='syn5730130'
csvs=afiles[grep("csv",afiles)]
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-02-17/testClusterParams.R'
ncats.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/ncatsSingleAgentScreens.R'
ctrp.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/ctrpSingleAgentScreens.R'
analysis.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/bin/singleDrugAnalysis.R'
for(csv in csvs){
  #first check to see if we're using CTRP or ncats, and original vs. rescored
  csv=paste('../2016-02-17',csv,sep='/')
  if(length(grep('ncats',csv))>0){
        
    if(length(grep('rescored',csv))>0){
      uf='syn5637634'
    }else{
      uf='syn5522627'
      
    }
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ncats.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE),
                          list(entity=uf,wasExecuted=FALSE)),
             activityName='drug AUC Clustering')
    
  }else if(length(grep('ctrp',csv))>0){
    if(length(grep('rescored',csv))>0){
      uf='syn5622708'
    }else{
      uf='syn5632189'
      
    }
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ctrp.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE),
                          list(entity=uf,wasExecuted=FALSE)),
             activityName='drug AUC Clustering')
    
  }else{
    sf=File(csv,parentId=cluster.dir)
    synStore(sf,used=list(list(url=this.script,wasExecuted=TRUE),
                          list(url=ncats.script,wasExecuted=TRUE),
                          list(url=ctrp.script,wasExecuted=TRUE),
                          list(url=analysis.script,wasExecuted=TRUE)),
             activityName='drug AUC Clustering Summary')
    
  }
}