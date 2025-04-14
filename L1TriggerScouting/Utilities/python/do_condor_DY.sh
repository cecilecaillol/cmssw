mkdir -p /eos/cms/store/group/cmst3/group/slowmuons/DYgen/KBMTFsim
mkdir -p submit_KBMTFsim_DY
python3 lxbatch_KBMTFsim.py --input_d=/eos/cms/store/group/cmst3/group/slowmuons/DYgen/RECOUSER --output_d=/eos/cms/store/group/cmst3/group/slowmuons/DYgen/KBMTFsim --process=DY
sh mysubmitKBMTFsim_DY.sh
