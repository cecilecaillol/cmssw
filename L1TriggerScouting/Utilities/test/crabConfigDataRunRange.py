from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'Scouting_2024G_run383804to384149'
#config.General.requestName = 'Scouting_2024G_run384151to384465'
#config.General.requestName = 'Scouting_2024G_run384467to384876'
#config.General.requestName = 'Scouting_2024G_run384877to385178'
#config.General.requestName = 'Scouting_2024G_run385179to385511'
config.General.requestName = 'Scouting_2024G_run385512to385813'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'kbmtFlatTableProducer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = "/L1Scouting/Run2024G-v1/L1SCOUT"
config.Data.inputDBS = "global"
config.Data.splitting = 'Automatic'
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json'
#config.Data.runRange = '383804-384149'
#config.Data.runRange = '384151-384465'
#config.Data.runRange = '384467-384876'
#config.Data.runRange = '384877-385178'
#config.Data.runRange = '385179-385511'
config.Data.runRange = '385512-385813'
config.Data.outLFNDirBase = '/store/group/cmst3/group/slowmuons/Mu8Skim/'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
