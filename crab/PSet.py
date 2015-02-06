import FWCore.ParameterSet.Config as cms

process = cms.Process('NoSplit')

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.output    = cms.OutputModule("PoolOutputModule",
    outputCommands  = cms.untracked.vstring("drop *", "keep recoTracks_*_*_*"),
    fileName        = cms.untracked.string('delphesTree.root'),
)
process.out       = cms.EndPath(process.output)

