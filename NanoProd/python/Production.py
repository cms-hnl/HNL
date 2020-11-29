# Produce NanoAOD customized for HNL studies.

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('sampleType', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Indicates the sample type: MC_18, Run2018ABC, ...")
options.register('fileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Text file with a list of root files to process.")
options.register('fileNamePrefix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Prefix to add to input file names.")
options.register('output', 'nanoTuple.root', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Output nanoAOD file.")
options.register('lumiFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "JSON file with lumi mask.")
options.register('eventList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "List of events to process.")
options.register('dumpPython', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Dump full config into stdout.")
options.register('storeExtInfo', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "If true, store extended information required for tau reconstruction")
options.parseArguments()


import HNL.NanoProd.sampleConfig as sampleConfig
isData = sampleConfig.isData(options.sampleType)
period = sampleConfig.getPeriod(options.sampleType)
process = sampleConfig.createProcess('NANO', options.sampleType)

process.options = cms.untracked.PSet()
process.options.wantSummary = cms.untracked.bool(False)

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

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))

# Set global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, sampleConfig.getGlobalTag(options.sampleType), '')

# Input and output definitions
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring())
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from HNL.NanoProd.readFileList import *
if len(options.fileList) > 0:
    readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
elif len(options.inputFiles) > 0:
    addFilesToList(process.source.fileNames, options.inputFiles, options.fileNamePrefix)

if options.maxEvents > 0:
    process.maxEvents.input = options.maxEvents

if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

if options.eventList != '':
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(',', options.eventList))

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(options.output),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Customisation of the process.

process.load('HNL.NanoProd.DiMuon_cff')
from HNL.NanoProd.DiMuon_cff import nanoAOD_customizeDisplacedDiMuon
nanoAOD_customizeDisplacedDiMuon(process)

from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC
process = nanoAOD_customizeMC(process)

process.schedule = cms.Schedule(
    process.nanoAOD_diDSAMuon_step,
    process.nanoAOD_patDSAMuon_step,
    process.nanoAOD_diMuon_step,
    process.endjob_step,
    process.NANOAODSIMoutput_step
)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# End of customisation functions

process.NANOAODSIMoutput.outputCommands = cms.untracked.vstring(
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

process.NANOAODSIMoutput.SelectEvents = cms.untracked.PSet(
    SelectEvents=cms.vstring(
        'nanoAOD_diDSAMuon_step',
        'nanoAOD_patDSAMuon_step',
        'nanoAOD_diMuon_step',
    )
)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

# Set report frequency
process.load('FWCore.MessageLogger.MessageLogger_cfi')
if process.maxEvents.input.value() > 10000 or process.maxEvents.input.value() < 0:
     process.MessageLogger.cerr.FwkReport.reportEvery = 1000
elif process.maxEvents.input.value() > 10:
     process.MessageLogger.cerr.FwkReport.reportEvery = process.maxEvents.input.value() // 10

if options.dumpPython:
    print process.dumpPython()