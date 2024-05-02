import os
import sys
import subprocess

process="Zjets"
doToy="false"
jobLabel="Job" #Job

path="/afs/cern.ch/user/d/dparedes/WorkCERN/Analysis_4topsSM/4topsSM/"
period=["all"]#"mc16a","mc16d","mc16e"]#,"all","mc16a","mc16d","mc16e"]

useWeight="true"
binning=["_Default","_plus"]
useBDT=["true","false"]

conf={}

#nJets
conf["_nJets_ge0"]="nJets>=0"
conf["_nJets_ge1"]="nJets>=1"

'''
conf["_Htlep_ge0_le200"]="HT_all-HT_jets<=200000"
conf["_Htlep_ge200_le400"]="(HT_all-HT_jets)>200000 && (HT_all-HT_jets)<=400000"
conf["_Htlep_ge400_le600"]="(HT_all-HT_jets)>400000 && (HT_all-HT_jets)<=600000"
conf["_Htlep_ge600_le800"]="(HT_all-HT_jets)>600000 && (HT_all-HT_jets)<=800000"
conf["_Htlep_ge0_le300"]="(HT_all-HT_jets)>0 && (HT_all-HT_jets)<=300000"
conf["_Htlep_ge300_le600"]="(HT_all-HT_jets)>300000 && (HT_all-HT_jets)<=600000"
conf["_Htlep_ge0_le400"]="(HT_all-HT_jets)>0 && (HT_all-HT_jets)<=400000"
conf["_Htlep_ge400_le800"]="(HT_all-HT_jets)>400000 && (HT_all-HT_jets)<=800000"
conf["_Htlep_ge800"]="(HT_all-HT_jets)>800000"
conf["_Htlep_ge600"]="(HT_all-HT_jets)>600000"
conf["_Htlep_ge400"]="(HT_all-HT_jets)>400000"

conf["_nJets_ge5"]="nJets>=5"
conf["_nJets_le4"]="nJets<=4"
conf["_nJets_ge4"]="nJets>=4"
conf["_nJets_le3"]="nJets<=3"
conf["_nJets_ge3"]="nJets>=3"
conf["_nJets_le2"]="nJets<=2"
'''

#Ht
'''
conf["_Ht_ge0_le200"]   ="HT_all>0 && HT_all<=200000"
conf["_Ht_ge200_le400"] ="HT_all>200000 && HT_all<=400000"
conf["_Ht_ge400_le600"] ="HT_all>400000 && HT_all<=600000"
conf["_Ht_ge600_le800"] ="HT_all>600000 && HT_all<=800000"
conf["_Ht_ge800_le1000"]="HT_all>800000 && HT_all<=1000000"
conf["_Ht_ge1000"]   ="HT_all>1000000"

conf["_Ht_ge0_le500"]="HT_all>0 && HT_all<=500000"
conf["_Ht_ge500_le1000"]="HT_all>500000 && HT_all<=1000000"
conf["_Ht_ge500"]   ="HT_all>500000"

conf["_Ht_ge0_le400"]="HT_all>0 && HT_all<=400000"
conf["_Ht_ge400_le800"]="HT_all>400000 && HT_all<=800000"
conf["_Ht_ge800"]   ="HT_all>800000"

conf["_Ht_ge0_le300"]="HT_all>0 && HT_all<=300000"
conf["_Ht_ge300_le600"]="HT_all>300000 && HT_all<=600000"
conf["_Ht_ge600"]   ="HT_all>600000"
conf["_Ht_ge600_le900"]="HT_all>600000 && HT_all<=900000"
conf["_Ht_ge900"]   ="HT_all>900000"
'''

'''
#nBJets
conf["_BJets_ge2"]="nBTags_MV2c10_77>=2"
conf["_BJets_le1"]="nBTags_MV2c10_77<=1"
conf["_BJets_ge3"]="nBTags_MV2c10_77>=3"
conf["_BJets_le2"]="nBTags_MV2c10_77<=2"

#MET
conf["_met_ge40_le80"] ="met_met>40000 && met_met<=80000"
conf["_met_le40"] ="met_met<40000"
conf["_met_ge80"] ="met_met>80000"
conf["_met_le80"] ="met_met<80000"
conf["_met_ge40"] ="met_met>40000"

#PV
conf["_pv_ge10"]="nPrimaryVtx>=10"
conf["_pv_ge20"]="nPrimaryVtx>=20"
conf["_pv_ge30"]="nPrimaryVtx>=30"

conf["_pv_ge0_le10"]="nPrimaryVtx>=0 && nPrimaryVtx<=10"
conf["_pv_ge10_le20"]="nPrimaryVtx>=10 && nPrimaryVtx<=20"
conf["_pv_ge20_le30"]="nPrimaryVtx>=20 && nPrimaryVtx<=30"
conf["_pv_ge0_le20"]="nPrimaryVtx>=0 && nPrimaryVtx<=20"

#HT_jets
conf["_Htjets_ge0_le200"]="HT_jets>0 && HT_jets<=200000"
conf["_Htjets_ge200_le400"]="HT_jets>200000 && HT_jets<=400000"
conf["_Htjets_ge400_le600"]="HT_jets>400000 && HT_jets<=600000"
conf["_Htjets_ge600_le800"]="HT_jets>600000 && HT_jets<=800000"
conf["_Htjets_ge0_le300"]="HT_jets>0 && HT_jets<=300000"
conf["_Htjets_ge300_le600"]="HT_jets>300000 && HT_jets<=600000"
conf["_Htjets_ge0_le400"]="HT_jets>0 && HT_jets<=400000"
conf["_Htjets_ge400_le800"]="HT_jets>400000 && HT_jets<=800000"
conf["_Htjets_ge800"]="HT_jets>800000"
conf["_Htjets_ge600"]="HT_jets>600000"
conf["_Htjets_ge400"]="HT_jets>400000"

#mu
conf["_mu_ge10"]="mu>=10"
conf["_mu_ge20"]="mu>=20"
conf["_mu_ge30"]="mu>=30"
conf["_mu_ge40"]="mu>=40"
conf["_mu_ge50"]="mu>=50"
conf["_mu_ge60"]="mu>=60"
conf["_mu_ge70"]="mu>=70"

conf["_mu_ge0_le30"]="mu>=0 && mu<=30"
conf["_mu_ge30_le40"]="mu>=30 && mu<=40"
conf["_mu_ge40_le55"]="mu>=40 && mu<=55"
conf["_mu_ge55"]="mu>=55"
conf["_mu_ge0_le20"]="mu>=0 && mu<=20"
conf["_mu_ge20_le40"]="mu>=20 && mu<=40"
conf["_mu_ge40_le60"]="mu>=40 && mu<=60"
'''

if not os.path.exists('scripts'):
    os.makedirs('scripts')

#for key, value in conf.iteritems():
#  print key, value

for per in period:

  if not os.path.exists('scripts/'+per):
    os.makedirs('scripts/'+per)
  if not os.path.exists('scripts/output/'+per):
    os.makedirs('scripts/output/'+per)
  if not os.path.exists('scripts/log/'+per):
    os.makedirs('scripts/log/'+per)
  if not os.path.exists('scripts/error/'+per):
    os.makedirs('scripts/error/'+per)

  for key, value in conf.iteritems():
    print key, value

  #for bin_local in binning:

    for BDT in useBDT:

      for bin_local in binning:
      #for per in period:
      
        queue='tomorrow'

        if 'mc16' in per:
          queue='tomorrow'


        script='''#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup "root 6.14.04-x86_64-slc6-gcc62-opt" 

cd {6}

pwd

root -l -q -b 'ChargeFlipBatchNoCF.C("{0}","{1}",{2},{3},"{4}","{5}","{7}",{8})'


'''.format(process,bin_local,useWeight,BDT,key,value,path,per,doToy)


        jobName=jobLabel+'_'+process+bin_local+key+"_"+per

        if BDT=="true":
            jobName+="_wBDT"

        if doToy=="true":
            jobName+="_Toy"
        else:
            jobName+="_Full"

        scriptFile=open('scripts/'+per+'/'+jobName+'.sh','w')
        scriptFile.write(script)
        scriptFile.close()

        os.chmod('scripts/'+per+'/'+jobName+'.sh',0o777)

        submitFile='''executable = {1}/{4}/{0}.sh
arguments  = $(ClusterID) $(ProcId)
output     = {1}/output/{4}/{0}.$(ClusterID).$(ProcId).out
error      = {1}/error/{4}/{0}.$(ClusterID).$(ProcId).err
log        = {1}/log/{4}/{0}.$(ClusterID).$(ProcId).log
request_disk ={3}

+JobFlavour = "{2}"
queue'''.format(jobName,path+'scripts/','workday','20000',per)

        scriptHTCondor=open('scripts/'+per+'/'+jobName+'_toSubmit.sh','w')
        scriptHTCondor.write(submitFile)
        scriptHTCondor.close()
        os.chmod('scripts/'+per+'/'+jobName+'_toSubmit.sh',0o777)

        bsubcommand='condor_submit '+path+'scripts/'+per+'/'+jobName+'_toSubmit.sh'

        print bsubcommand

        
        subprocess.call(bsubcommand, shell=True)

