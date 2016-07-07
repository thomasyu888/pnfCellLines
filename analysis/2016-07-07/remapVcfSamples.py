'''
Grab VCF file and re-header with proper sample names
'''


import os,re,sys

##first get the vcf file
os.system('synapse get syn6091968')


##now just get the samples
os.system('bcftools query -l H_Blakeley_NF1_WGHum-SeqWholeExome_160219_1.BEDsuperset.BestPractices.HardCutoffs.vcf>sampnames.txt')

##now read in the sample names
samps=[a.strip('\n') for a in open('sampnames.txt').readlines()]

import synapseclient
syn = synapseclient.Synapse()
syn.login()

##now get sample names to cell line names
tq="SELECT 'Sample Name','Exome-Seq Identifier' FROM syn5014742 where \"Exome-Seq Identifier\" is not NULL"


tq=syn.tableQuery(tq)

tab=[a.strip('\n').split(',')[2:4] for a in open(tq.filepath).readlines()[1:]]

newsamps=[]
for s in samps:
    tup=[a[0].strip() for a in tab if a[1].strip('"')==s]
    print tup
    newsamps.append(tup[0].strip('"')+'\n')

open('newHeader.txt','w').writelines(newsamps)

##nowre-header
newfile= 'H_Blakeley_NF1_WGHum-SeqWholeExome_160219_1.BEDsuperset.BestPractices.HardCutoffs_WithSamples.vcf'
os.system('bcftools reheader -s newHeader.txt H_Blakeley_NF1_WGHum-SeqWholeExome_160219_1.BEDsuperset.BestPractices.HardCutoffs.vcf>'+newfile)

##gzip
os.system('bgzip '+newfile)
##now index
os.system('bcftools index -t '+newfile+'.gz')
