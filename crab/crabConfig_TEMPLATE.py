from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName  = 'SAMPLENAME_FIXME'
config.General.workArea     = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName   = 'PrivateMC'
config.JobType.psetName     = 'FOLDER_FIXME/PSet.py'
config.JobType.scriptExe    = 'FOLDER_FIXME/delphes_job.sh'
config.JobType.inputFiles   = [
        'FOLDER_FIXME/FrameworkJobReport.xml', 
        'FOLDER_FIXME/Delphes.tgz', 
        'FOLDER_FIXME/fastjet-3.1.0.tar.gz', 
        'FOLDER_FIXME/pythia8201.tgz', 
        'FOLDER_FIXME/inputlist.txt', 
        'FOLDER_FIXME/pulist.txt', 
        'FOLDER_FIXME/the_fjcontrib.tgz', 
#        'FOLDER_FIXME/CMS_Phase_I_140PileUp_Tracker2p5.tcl',
#        'FOLDER_FIXME/PROXY_FIXME'
    ]
#config.JobType.outputFiles  = ['output.txt']
config.JobType.outputFiles  = ['delphesTree.root']
#config.JobType.outputFiles  = ['output.txt', 'delphesTree.root']


config.section_("Data")
#config.Data.inputDataset    = '/SingleMu/Run2012B-13Jul2012-v1/AOD'
#config.Data.inputDBS        = 'global'
config.Data.primaryDataset = 'ttbar'
config.Data.splitting      = 'EventBased'
#config.Data.splitting       = 'LumiBased'
config.Data.unitsPerJob = 1 
NJOBS =   NJOBS_FIXME                # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.lumiMask        = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange        = '193093-193999' # '193093-194075'
config.Data.outLFN          = '/store/user/mmozer/Delphes/'
config.Data.publication     = False
# lsconfig.Data.trasferoutput   = False

config.section_("Site")
config.Site.storageSite     = 'T2_DE_DESY'
#config.Site.whitelist       = ['T2_CH_CERN','T2_RU_JINR']
#config.Site.storageSite     = 'T3_IT_MIB'
#config.Site.storageSite     = 'T1_DE_KIT'
