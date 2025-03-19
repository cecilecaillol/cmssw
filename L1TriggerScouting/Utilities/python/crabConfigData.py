from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'Scouting_2024G'
config.General.requestName = 'Scouting_2024H'
#config.General.requestName = 'Scouting_2024I'
#config.General.requestName = 'Scouting_2024J'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'kbmtFlatTableProducer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

#config.Data.inputDataset = "/L1Scouting/Run2024G-v1/L1SCOUT"
config.Data.inputDataset = "/L1Scouting/Run2024H-v1/L1SCOUT"
#config.Data.inputDataset = "/L1Scouting/Run2024I-v1/L1SCOUT"
#config.Data.inputDataset = "/L1Scouting/Run2024J-v1/L1SCOUT"
config.Data.inputDBS = "global"
config.Data.splitting = 'Automatic'
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json'
config.Data.outLFNDirBase = '/store/group/cmst3/group/slowmuons/Mu8Skim/'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
