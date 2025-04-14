#mkdir -p /eos/cms/store/group/cmst3/group/slowmuons/DYgen/KBMTFsim
#mkdir -p submit_KBMTFsim_DY
#python3 lxbatch_KBMTFsim.py --input_d=/eos/cms/store/group/cmst3/group/slowmuons/DYgen/RECOUSER --output_d=/eos/cms/store/group/cmst3/group/slowmuons/DYgen/KBMTFsim --process=DY
#sh mysubmitKBMTFsim_DY.sh

#!/bin/bash

masses=("2600" "3000" "3500" "4000" "4500" "5000" "5500" "6000")

for mass in "${masses[@]}"; do
   mkdir -p /eos/cms/store/group/cmst3/group/slowmuons/HSCPtauPrime_${mass}/KBMTFsim
   mkdir -p submit_KBMTFsim_HSCPtauPrime_${mass}
   python3 lxbatch_KBMTFsim.py --input_d=/eos/cms/store/group/cmst3/group/slowmuons/HSCPtauPrime_${mass}/RECOUSER --output_d=/eos/cms/store/group/cmst3/group/slowmuons/HSCPtauPrime_${mass}/KBMTFsim --process=HSCPtauPrime_${mass}
   sh mysubmitKBMTFsim_HSCPtauPrime_${mass}.sh
done

