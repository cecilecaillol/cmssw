import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_d")
parser.add_argument("--output_d")
parser.add_argument("--process")
input_d = parser.parse_args().input_d
output_d = parser.parse_args().output_d
process = parser.parse_args().process

list_f = os.listdir(input_d)#.listdir()

submit_command = open("mysubmitKBMTFsim_"+process+".sh","w")
sub_file = open("submission_KBMTFsim_"+process+".sub","w")
sub_file.write("executable            = submit_KBMTFsim_"+process+"/hello_world_$(ProcId).sh \n")
sub_file.write("output                = output_KBMTFsim_"+process+"/hello.$(ProcId).out \n")
sub_file.write("error                 = error_KBMTFsim_"+process+"/hello.$(ProcId).err \n")
sub_file.write("log                   = log_KBMTFsim_"+process+"/hello.$(ProcId).log \n")
sub_file.write('+JobFlavour = "workday" \n')
sub_file.write('+AccountingGroup = "group_u_CMST3.all" \n')
sub_file.write("queue "+str(len(list_f))+"\n")

submit_command.write("mkdir submit_KBMTFsim_"+process+" \n")
submit_command.write("mkdir output_KBMTFsim_"+process+" \n")
submit_command.write("mkdir error_KBMTFsim_"+process+" \n")
submit_command.write("mkdir log_KBMTFsim_"+process+" \n")

for i in range(0,len(list_f)):
  exe_file = open("submit_KBMTFsim_"+process+"/hello_world_"+str(i)+".sh","w")

#!/bin/bash
  exe_file.write('#!/bin/bash \n')
  exe_file.write('source /cvmfs/cms.cern.ch/cmsset_default.sh \n')
  exe_file.write('cd /afs/cern.ch/work/c/ccaillol/KBMTFemulation/CMSSW_14_0_12/src/L1TriggerScouting/Utilities/python \n')
  exe_file.write('eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers \n')
  exe_file.write('echo $CMSSW_BASE "is the CMSSW we created on the local worker node" \n')
  exe_file.write('cd submit_KBMTFsim_'+process +' \n')
  exe_file.write("cmsRun run_KBMTFsim_"+str(i)+".py \n ")
  exe_file.close()
  submit_command.write("python3 replace.py --skeleton skeleton_kbmtFlatTableProducerSimulation_cfg.py --output_py submit_KBMTFsim_"+process+"/run_KBMTFsim_"+str(i)+".py --input_root "+input_d+"/"+list_f[i]+" --output_root "+output_d+"/KBMTFsim_"+list_f[i]+" \n\n") 

submit_command.write("condor_submit submission_KBMTFsim_"+process+".sub")
  
submit_command.close()
