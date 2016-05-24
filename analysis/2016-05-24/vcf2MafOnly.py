'''
This is the primary script that can be used to process the VCFs after
receiving data from CIDR
'''

##VCF2MAF installation:
vcf2maf='perl ../../../vcf2maf-master/vcf2maf.pl'
##REFERENCE FASTA installation
##UPDATE WITH new
reffasta='../../lib/ucsc.hg19.fasta'
##OUTPUT FROM GATK:



import synapseclient,re,os
syn = synapseclient.Synapse()
syn.login()

vcf_file=syn.get('syn5555584').path


##get all the VCF annotations so that we can process the merged file
wgs_vcf='syn6086903'
#get all files
query_res=syn.query("select id, name from entity where entity.parentId=='"+wgs_vcf+"'")
##now get metadata for all files
syn_files=dict([(res['entity.id'],res['entity.name']) for res in query_res['results'] if 'idx' not in res['entity.name']])
syn_ids=syn_files.keys()

##if this is what we're going to use as control, try that.
control=[k for k in syn_files.keys() if 'NA12891' in syn_files[k]]
control_vcf=syn.get(str(control[0])).path

##now get metadata only for these files to determine which synapse ids belong to
##a specific patient
#annotes=[syn.getAnnotations(a) for a in syn_ids]
annotes=syn.tableQuery('SELECT "Sample Name","Exome-Seq Identifier","Sample Origin","Sample Genotype" FROM syn5014742 where "Exome-Seq Identifier" is not NULL').asDataFrame()

for row in annotes.iterrow():
    ##get all the values
    exomeid=row[1][1]
    if 'NA12891' in exomeid:
        continue

    orig=row[1][2]
    gt=row[1][3]
    sampname=row[1][0]
    nospace=re.sub('(|)','_',re.sub(' ','_',sampname))
    print nospace
    samp_id=[k for k in syn_files.keys() if exomeid in syn_files[k]]
    samp_vcf=syn.get(str(samp_id[0])).path


    cmfile='sample_'+nospace+'_commands.sh'
        ##first create vcf file with only two samples
#    bcftoolscmd="bcftools view %s -s %s -U"%(samp_vcf,sampname)
#    outvcf="sample_%s_normOnly_%s.vcf"%(nospace,samp_id)
    outmaf=re.sub('.vcf','.maf',os.path.basename(samp_vcf))

                    #create new annotation string
        #activity string?
#    annotationstr="'{\"dataType\":\"WGS\",\"tissueType\":\"PBMC\",\"patientId\":\""+bloodfile_annotations['patientID'][0]+"\"}'"
#        annotationstr=annotationstr+",\"tissueID\":\""+tumfile_annotations['tissueID'][0]+"\"}'"



#	if bloodfile=='':
#	    usedstr="'{"+t+"}'"
#	else:
#	    usedstr="'{"+t+","+blood[0]+"}'"
        ##now store paired VCF
        #synapse_upload_vcf="synapse store "+outvcf+".gz --parentId=syn5522791 --annotations "+annotationstr#+' --used '+usedstr

        #patsh.write(synapse_upload_vcf+'\n\n')

        #patsh.write('bgzip -d '+outvcf+'.gz\n\n')
    if os.path.exists(outmaf+'.gz'):
        print outmaf+'.gz is already made, not re-creating'
        continue

    patsh=open(cmfile,'w')
    vcf2maf_cmd=vcf2maf+" --input-vcf %s --vcf-normal-id %s"%(samp_vcf,sampname)
    vcf2maf_cmd+=" --output-maf %s --vep-forks 16 --species homo_sapiens --ref-fasta %s"%(outmaf,reffasta)
    patsh.write(bcftoolscmd+'>'+outvcf+'\n') #bgzip '+outvcf+'\n')

    patsh.write(vcf2maf_cmd+'\ngzip '+outmaf+'\n')
        #then these file should be uploaded to synapse

    synapse_upload_maf="synapse store "+outmaf+".gz --parentId=syn6117796"
#--annotations "+annotationstr#+' --used '+usedstr
    patsh.write(synapse_upload_maf+'\n')

        #patsh.write(bcftoolscmd+'>'+outvcf+'\n'+vcf2maf_cmd+'\n')
        #os.system(cmd)
    patsh.close()
    #os.system('sh '+cmfile+' &')

        #examples for patient 1
