import os
import sys
import subprocess




outputDir="ValidationPlots_212560"


jobLabel="Job" #Job

path="/afs/cern.ch/user/d/dparedes/WorkCERN/Analysis_4topsSM/4topsSM/"



#if process=="Zjets" or process=="ttbar":  #do normal splitting just when checking truth rates. Otherwise period=process. 
#    period=["mc","mc16a","mc16d","mc16e"]


useBDT=["true"]#,"false"]
variables=["pt1","pt2","eta1","met","Ht","Htlep","Htjets","nJets","BJets","mu","pv","Mll"]
testHt=["false"]#["true","false"]

if not os.path.exists('scriptsValidation'):
    os.makedirs('scriptsValidation')

#for key, value in conf.iteritems():
#  print key, value


process=["Data","Data2015","Data2016","Data2017","Data2018","mc","mc16a","mc16d","mc16e"]  #ttbar,Zjets,mc,Data


for proc in process:
   
    period=[proc]

    for per in period:

      if not os.path.exists('scriptsValidation/'+per):
	os.makedirs('scriptsValidation/'+per)
      if not os.path.exists('scriptsValidation/output/'+per):
	os.makedirs('scriptsValidation/output/'+per)
      if not os.path.exists('scriptsValidation/log/'+per):
	os.makedirs('scriptsValidation/log/'+per)
      if not os.path.exists('scriptsValidation/error/'+per):
	os.makedirs('scriptsValidation/error/'+per)

      queue='longlunch'

      #if 'mc16' in per:
      #  queue='workday'


      for var in variables:
	print var

	for BDT in useBDT:

	  for test in testHt:

            if 'false' in BDT:
              inputDir="Output212560_Weighted_Full"         
            else:
              inputDir="Output212560_Weighted_wBDT_Full"
            localOutput=outputDir

            if 'true' in test:
              localOutput+="_testH"


    #"ttbar", "Output_Weighted_Full", "mc16d", "v_Ht", "ValidationTest"))
            script='''#!/bin/bash

    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
    lsetup "root 6.14.04-x86_64-slc6-gcc62-opt" 

    cd {5}

    pwd

    root -l -q -b 'MethodABatchFull.C("{0}","{1}","{2}","{3}","{4}",{6})'


    '''.format(proc,inputDir,per,var,localOutput,path,test)


            jobName=jobLabel+'_'+proc+"_"+var+"_"+per

            if BDT=="true":
        	jobName+="_wBDT"

            if test=="true":
        	jobName+="_testHt"

            scriptFile=open('scriptsValidation/'+per+'/'+jobName+'.sh','w')
            scriptFile.write(script)
            scriptFile.close()

            os.chmod('scriptsValidation/'+per+'/'+jobName+'.sh',0o777)

            submitFile='''executable = {1}/{4}/{0}.sh
arguments  = $(ClusterID) $(ProcId)
output     = {1}/output/{4}/{0}.$(ClusterID).$(ProcId).out
error      = {1}/error/{4}/{0}.$(ClusterID).$(ProcId).err
log        = {1}/log/{4}/{0}.$(ClusterID).$(ProcId).log
request_disk ={3}

+JobFlavour = "{2}"
queue'''.format(jobName,path+'/scriptsValidation/','workday','20000',per)

            scriptHTCondor=open('scriptsValidation/'+per+'/'+jobName+'_toSubmit.sh','w')
            scriptHTCondor.write(submitFile)
            scriptHTCondor.close()
            os.chmod('scriptsValidation/'+per+'/'+jobName+'_toSubmit.sh',0o777)

            bsubcommand='condor_submit '+path+'scriptsValidation/'+per+'/'+jobName+'_toSubmit.sh'

            print bsubcommand


            subprocess.call(bsubcommand, shell=True)

