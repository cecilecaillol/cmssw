#mkdir -p /eos/cms/store/group/cmst3/group/slowmuons/DYgen/KBMTFsim
#mkdir -p submit_KBMTFsim_DY
#python3 lxbatch_KBMTFsim.py --input_d=/eos/cms/store/group/cmst3/group/slowmuons/DYgen/RECOUSER --output_d=/eos/cms/store/group/cmst3/group/slowmuons/DYgen/KBMTFsim --process=DY
#sh mysubmitKBMTFsim_DY.sh

#!/bin/bash


for pair in "3000 1500" "4000 1500" "4000 2000" "5000 2000" "5000 2500" "6000 2500" "6000 3000" "7000 2500" "7000 3000" "7000 3500" "8000 3000" "8000 3500" "8000 4000" "9000 3500" "9000 4000" "9000 4500" "10000 4000" "10000 4500" "10000 5000" "11000 4500" "11000 5000" "11000 5500" "12000 4500" "12000 5000" "12000 5500" "12000 6000" ; do
   set -- $pair
   mZ=$1
   mtau=$2
   mkdir -p /eos/cms/store/group/cmst3/group/slowmuons/ZPrimeTo2TauPrime_${mZ}_${mtau}/KBMTFsim
   mkdir -p submit_KBMTFsim_ZPrimeTo2TauPrime_${mZ}_${mtau}
   python3 lxbatch_KBMTFsim.py --input_d=/eos/cms/store/group/cmst3/group/slowmuons/ZPrimeTo2TauPrime_${mZ}_${mtau}/RECOUSER --output_d=/eos/cms/store/group/cmst3/group/slowmuons/ZPrimeTo2TauPrime_${mZ}_${mtau}/KBMTFsim --process=ZPrimeTo2TauPrime_${mZ}_${mtau}
   sh mysubmitKBMTFsim_ZPrimeTo2TauPrime_${mZ}_${mtau}.sh
done

