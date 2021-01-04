# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: nano_v7_data_cff -s NANO --eventcontent NANOAOD --datatier NANOAOD --no_exec --conditions 102X_dataRun2_v13 --era Run2_2018,run2_nanoAOD_102Xv1
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_102Xv1)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HNL.NanoProd.DiMuon_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/90000/6FEE6A9B-E03C-A14F-AD01-914C76D0A721.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('nano_v7_data_cff nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('nano_v7_data_cff_NANO.root'),
    outputCommands = process.NANOAODEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v13', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# # customisation of the process.

from HNL.NanoProd.DiMuon_cff import nanoAOD_customizeDisplacedDiMuon
nanoAOD_customizeDisplacedDiMuon(process)

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeData(process)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_diDSAMuon_step, process.nanoAOD_patDSAMuon_step, process.nanoAOD_diMuon_step, process.endjob_step, process.NANOAODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# End of customisation functions

process.NANOAODoutput.outputCommands = cms.untracked.vstring(
    'drop *', 
    'keep nanoaodFlatTable_*Table_*_*', 
    'drop nanoaodFlatTable_fatJet*Table_*_*',
    'drop nanoaodFlatTable_subJet*Table_*_*',
    'drop nanoaodFlatTable_subjet*Table_*_*',
    'drop nanoaodFlatTable_softActivity*Table_*_*',
    'drop nanoaodFlatTable_farPhotonTable_*_*',
    'drop nanoaodFlatTable_HTXSCategoryTable_*_*',
    'drop nanoaodFlatTable_corrT1METJetTable_*_*',
    'drop nanoaodFlatTable_sa*Table_*_*',
    'drop nanoaodFlatTable_svCandidateTable_*_*',
    'drop nanoaodFlatTable_rivet*Table_*_*',
    # it does not seem possible to skim the trigger information saved with existing tools
    'keep edmTriggerResults_*_*_*',  
    'keep String_*_genModel_*', 
    'keep nanoaodMergeableCounterTable_*Table_*_*', 
    'keep *MergeableCounterTable_*_*_*', 
    'keep nanoaodUniqueString_nanoMetadata_*_*'
)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
    SelectEvents=cms.vstring(
        'nanoAOD_diDSAMuon_step',
        'nanoAOD_patDSAMuon_step',
        'nanoAOD_diMuon_step',
    )
)

# End of customisation functions

# Customisation from command line
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(20)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
