##upload files
library(synapseClient)
synapseLogin()
synid='syn5731320'

afiles=list.files("../2016-01-14")
this.script='https://raw.githubusercontent.com/Sage-Bionetworks/NTAP/master/pnfCellLines/analysis/2016-01-14/testSleuth.R'

ufiles=afiles[-grep('.R',afiles)]
for(f in ufiles){
    synStore(File(paste('../2016-01-14',f,sep='/'),parentId=synid),used=list(list(url=this.script,executed=TRUE)))
  
  
}