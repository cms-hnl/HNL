import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, Var, P3Vars
from HNL.NanoProd.common_cff import ufloat, uint, ubool

from HNL.NanoProd.DiMuon_cff import *

pfCandidateTable = cms.EDProducer('SimpleCandidateFlatTableProducer',
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string(""),
    name = cms.string("pfCand"),
    doc = cms.string("PF candidates"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        P3Vars,
        charge = Var("charge", int, doc="electric charge"),
        pdgId = Var("pdgId", int, doc="pdg ID"),
        dxy = Var('dxy', float, precision=10, doc='dxy'),
        dz = Var('dz', float, precision=10, doc='dz'),
        vertex_x = Var('vertex().x', float, precision=10, doc='dxy'),
        vertex_y = Var('vertex().y', float, precision=10, doc='dxy'),
        vertex_z = Var('vertex().z', float, precision=10, doc='dxy'),
    )
)

def nanoAOD_customizeDisplacedTau(process):
    process.displacedDiMuonSequence = cms.Sequence(diDSAMuonSequence*patDSAMuonSequence*diMuonSequence)
    # process.nanoAOD_step.insert(1000, process.displacedDiMuonSequence)
    process.muonSequence.insert(1000, process.displacedDiMuonSequence)
    process.nanoAOD_diDSAMuon_step = cms.Path(process.nanoSequenceMC + countDiDSAMuon)
    process.nanoAOD_patDSAMuon_step = cms.Path(process.nanoSequenceMC + countPatDSAMuon)
    process.nanoAOD_diMuon_step = cms.Path(process.nanoSequenceMC + countDiMuon)

    process.finalMuons.cut = "pt > 3"

    # Add additional muon time variables
    process.muonTable.variables.rpcTimeInOut = Var("rpcTime().timeAtIpInOut", float, doc="RPC time in out")
    process.muonTable.variables.timeInOut = Var("time().timeAtIpInOut", float, doc="time in out")
    process.muonTable.variables.rpcTimeInOutErr = Var("rpcTime().timeAtIpInOutErr", float, doc="RPC time error in out")
    process.muonTable.variables.timeInOutErr = Var("time().timeAtIpInOutErr", float, doc="time error in out")
    process.muonTable.variables.rpcTimeNdof = Var("rpcTime().nDof", float, doc="RPC time ndof")
    process.muonTable.variables.timeNdof = Var("time().nDof", float, doc="time ndof")



    process.isoTrackTable.externalVariables = cms.PSet()
    process.isoTrackTable.src = cms.InputTag("isolatedTracks")

    process.lostTrackTable = process.isoTrackTable.clone(
        src = cms.InputTag("lostTracks", ""),
        name = cms.string("LostTrack")
    )
    for v in [ "isHighPurityTrack", "isPFcand", "pfRelIso03_all", "pfRelIso03_chg" ]:
        delattr(process.lostTrackTable.variables, v)

    process.lostTrackEleTable = process.lostTrackTable.clone(
        src = cms.InputTag("lostTracks", "eleTracks"),
        name = cms.string("LostTrackEle")
    )

    process.nanoSequenceMC.insert(-1, process.lostTrackTable)
    process.nanoSequenceMC.insert(-1, process.lostTrackEleTable)
    process.nanoSequenceMC.insert(-1, process.pfCandidateTable)

    process.genParticleTable.variables.vertex_x = Var('vertex().x', float, precision=10, doc='dxy')
    process.genParticleTable.variables.vertex_y = Var('vertex().y', float, precision=10, doc='dxy')
    process.genParticleTable.variables.vertex_z = Var('vertex().z', float, precision=10, doc='dxy')
    process.genParticleTable.variables.vertex_rho = Var('vertex().rho', float, precision=10, doc='dxy')
    process.genParticleTable.variables.vertex_r = Var('vertex().r', float, precision=10, doc='dxy')


    # Skim jet variables
    process.jetTable.externalVariables = cms.PSet()
    jet_vars_rm = ['btagCMVA', 'btagDeepB', 'btagCSVV2', 'qgl',  'chHEF', 'neHEF', 'chEmEF', 'neEmEF', 'muEF',
                   'chFPV0EF', 'chFPV1EF', 'chFPV2EF', 'chFPV3EF']
    for v in jet_vars_rm:
        delattr(process.jetTable.variables, v)

    # remove HTXS table
    process.particleLevelSequence.remove(process.rivetProducerHTXS)
    process.particleLevelTables.remove(process.HTXSCategoryTable)

    # remove fat jets
    process.jetTables.remove(process.fatJetTable)
    process.nanoSequenceCommon.remove(process.jetLepSequence)
    process.jetMC.remove(process.fatJetMCTable)
    for c in [ 'jetCorrFactorsAK8', 'updatedJetsAK8', 'tightJetIdAK8', 'tightJetIdLepVetoAK8',
              'updatedJetsAK8WithUserData', 'finalJetsAK8' ]:
        process.jetSequence.remove(getattr(process, c))
    import PhysicsTools.NanoAOD.nano_cff as nano_cff
    def fake_add(process, *args, **kwargs):
        return process
    nano_cff.nanoAOD_addDeepInfoAK8 = fake_add

    return process
