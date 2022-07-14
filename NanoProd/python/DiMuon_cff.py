import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, Var, P3Vars
from HNL.NanoProd.common_cff import ufloat, uint, ubool

selectedDSAMuons = cms.EDFilter(
    'TrackSelector',
    src = cms.InputTag('displacedStandAloneMuons'),
    cut = cms.string('pt > 5. && abs(eta) < 2.4 && numberOfValidHits > 15 && ptError/pt < 1. && chi2/ndof < 2.5')
)

vetoMuons = cms.EDFilter(
    'PATMuonRefSelector',
    src = cms.InputTag('slimmedMuons'),
    cut = cms.string('pt>24 && abs(eta) < 2.4 && isMediumMuon && dB<0.02')
)


diDSAMuon = cms.EDProducer(
    'DiTrackBuilder',
    src = cms.InputTag('selectedDSAMuons'),
    srcVeto = cms.InputTag('vetoMuons'),
    lep1Selection = cms.string('1'),
    lep2Selection = cms.string('1'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('userFloat("sv_ndof") > 0')
)

patDSAMuon = cms.EDProducer(
    'MuTrackBuilder',
    src1 = cms.InputTag('finalMuons'),
    src2 = cms.InputTag('selectedDSAMuons'),
    srcVeto = cms.InputTag('vetoMuons'),
    lep1Selection = cms.string('pt > 3. && isGlobalMuon && dB > 0.01 && isLooseMuon && abs(eta) < 2.4 && segmentCompatibility > 0.451 && combinedQuality().trkKink < 20'),
    lep2Selection = cms.string('1'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('userFloat("sv_ndof") > 0')
)

diMuon = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('finalMuons'),
    srcVeto = cms.InputTag('vetoMuons'),
    lep1Selection = cms.string('pt > 3. && isGlobalMuon && dB > 0.01 && isLooseMuon && abs(eta) < 2.4 && segmentCompatibility > 0.451 && combinedQuality().trkKink < 20'),
    lep2Selection = cms.string('pt > 3. && isGlobalMuon && dB > 0.01 && isLooseMuon && abs(eta) < 2.4 && segmentCompatibility > 0.451 && combinedQuality().trkKink < 20'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('userFloat("sv_ndof") > 0')
)

diSTA = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('finalMuons'),
    srcVeto = cms.InputTag('vetoMuons'),
    lep1Selection = cms.string('pt > 5. && isStandAloneMuon && abs(eta) < 2.4 && numberOfValidHits > 12 && bestTrack.ptError/pt < 1. && bestTrack.chi2/bestTrack.ndof < 2.5'),
    lep2Selection = cms.string('pt > 5. && isStandAloneMuon && abs(eta) < 2.4 && numberOfValidHits > 12 && bestTrack.ptError/pt < 1. && bestTrack.chi2/bestTrack.ndof < 2.5'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('userFloat("sv_ndof") > 0')
)

patSTA = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('finalMuons'),
    srcVeto = cms.InputTag('vetoMuons'),
    lep1Selection = cms.string('pt > 5. && isStandAloneMuon && abs(eta) < 2.4 && numberOfValidHits > 12 && bestTrack.ptError/pt < 1. && bestTrack.chi2/bestTrack.ndof < 2.5'),
    lep2Selection = cms.string('pt > 3. && isGlobalMuon && dB > 0.01 && isLooseMuon && abs(eta) < 2.4 && segmentCompatibility > 0.451 && combinedQuality().trkKink < 20'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('userFloat("sv_ndof") > 0')
)


dsaTable = cms.EDProducer(
    'SimpleTrackFlatTableProducer',
    src = cms.InputTag("selectedDSAMuons"),
    cut = cms.string("1"), # if we place a cut here, the indexing will be wrong
    name = cms.string("DSAMuon"),
    doc = cms.string("Displaced standalone muon tracks variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        P3Vars,
        charge = Var("charge", int, doc="electric charge"),
        n_valid_hits = Var('numberOfValidHits', int, doc='valid hits'),
        n_lost_hits = Var('numberOfLostHits', int, doc='lost hits'),
        n_muon_stations = Var('hitPattern().muonStationsWithValidHits', int, doc='muon stations with valid hits'),
        n_dt_stations = Var('hitPattern().dtStationsWithValidHits', int, doc='DT stations with valid hits'),
        n_dt_hits = Var('hitPattern().numberOfValidMuonDTHits', int, doc='valid DT hits'),
        n_csc_stations = Var('hitPattern().cscStationsWithValidHits', int, doc='CSC stations with valid hits'),
        n_csc_hits = Var('hitPattern().numberOfValidMuonCSCHits', int, doc='valid CSC hits'),
        n_rpc_stations = Var('hitPattern().rpcStationsWithValidHits', int, doc='RPC stations with valid hits'),
        n_rpc_hits = Var('hitPattern().numberOfValidMuonRPCHits', int, doc='valid RPC hits'),
        chi2 = Var('chi2', float, precision=10, doc='track chi2'),
        ndof = Var('ndof', float,precision=10,  doc='track ndof'),
        dxy = Var('dxy', float, precision=10, doc='dxy'),
        dz = Var('dz', float, precision=10, doc='dz'),
        pt_error = Var('ptError', float, precision=10, doc='pt error'),
        theta_error = Var('thetaError', float, precision=8, doc='theta error'),
        phi_error = Var('phiError', float, precision=8, doc='phi error'),
    )
)

dsaIsoTable = cms.EDProducer(
    'TrackIsoTableProducer',
    name = cms.string("DSAMuon")
)

diDSAMuonTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("diDSAMuon"),
    cut = cms.string(""),
    name = cms.string("DiDSAMuon"),
    doc = cms.string("DiDSAMuon Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),

        # fit and vtx info
        chi2 = ufloat('sv_chi2', precision=8),
        ndof = ufloat('sv_ndof', precision=8),
        svprob = ufloat('sv_prob', precision=10),
        # l_xy = ufloat('l_xy'),
        # l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x', precision=10),
        vtx_y = ufloat('vtx_y', precision=10),
        vtx_z = ufloat('vtx_z', precision=10),
        vtx_ex = ufloat('vtx_ex', precision=10), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey', precision=10),
        vtx_ez = ufloat('vtx_ez', precision=10),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass', precision=8),
        fit_massErr = ufloat('fitted_massErr', precision=8),
        fit_pt = ufloat('fitted_pt', precision=8),
    )
)

patDSAMuonTable = diDSAMuonTable.clone(
    src='patDSAMuon',
    name='PatDSAMuon',
    doc='PatDSAMuon Variable'
)

diMuonTable = diDSAMuonTable.clone(
    src='diMuon',
    name='DiMuon',
    doc='DiMuon Variable'
)

diSTATable = diDSAMuonTable.clone(
    src='diSTA',
    name='DiSTA',
    doc='DiSTA Variable'
)

patSTATable = diDSAMuonTable.clone(
    src='patSTA',
    name='PatSTA',
    doc='PatSTA Variable'
)


countDiDSAMuon = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("diDSAMuon")
)

countPatDSAMuon = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("patDSAMuon")
)

countDiMuon = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("diMuon")
)

countDiSTA = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("diSTA")
)

countPatSTA = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("patSTA")
)

diDSAMuonTables = cms.Sequence(dsaTable*diDSAMuonTable)
diDSAMuonSequence = cms.Sequence(selectedDSAMuons*vetoMuons*diDSAMuon*dsaTable*dsaIsoTable*diDSAMuonTable)
patDSAMuonSequence = cms.Sequence(patDSAMuon*patDSAMuonTable)
diMuonSequence = cms.Sequence(diMuon*diMuonTable)
diSTASequence = cms.Sequence(diSTA*diSTATable)
patSTASequence = cms.Sequence(patSTA*patSTATable)


isomu24 = cms.EDFilter('TriggerResultsFilter',
    hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
    l1tResults = cms.InputTag(''),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    throw = cms.bool(True),
    triggerConditions = cms.vstring('HLT_IsoMu24_v*'),
    # mightGet = cms.optional.untracked.vstring
)

def nanoAOD_customizeDisplacedDiMuon(process, is_mc=False):
    process.displacedDiMuonSequence = cms.Sequence(diDSAMuonSequence*patDSAMuonSequence*diMuonSequence*diSTASequence*patSTASequence)
    # process.nanoAOD_step.insert(1000, process.displacedDiMuonSequence)
    process.muonSequence.insert(1000, process.displacedDiMuonSequence)
    nano_seq = process.nanoSequenceMC if is_mc else process.nanoSequence
    process.nanoAOD_diDSAMuon_step = cms.Path(nano_seq + countDiDSAMuon)
    process.nanoAOD_patDSAMuon_step = cms.Path(nano_seq + isomu24 + countPatDSAMuon)
    process.nanoAOD_diMuon_step = cms.Path(nano_seq + isomu24 + countDiMuon)
    process.nanoAOD_diSTA_step = cms.Path(nano_seq + isomu24 + countDiSTA)
    process.nanoAOD_patSTA_step = cms.Path(nano_seq + isomu24 + countPatSTA)

    process.finalMuons.cut = "pt > 3"

    # Add additional muon time variables
    process.muonTable.variables.rpcTimeInOut = Var("rpcTime().timeAtIpInOut", float, doc="RPC time in out")
    process.muonTable.variables.timeInOut = Var("time().timeAtIpInOut", float, doc="time in out")
    process.muonTable.variables.rpcTimeInOutErr = Var("rpcTime().timeAtIpInOutErr", float, doc="RPC time error in out")
    process.muonTable.variables.timeInOutErr = Var("time().timeAtIpInOutErr", float, doc="time error in out")
    process.muonTable.variables.rpcTimeNdof = Var("rpcTime().nDof", float, doc="RPC time ndof")
    process.muonTable.variables.timeNdof = Var("time().nDof", float, doc="time ndof")
    process.muonTable.variables.nValidHits = Var("numberOfValidHits", float, doc="n valid hits")
    process.muonTable.variables.chi2ndof = Var("bestTrack.chi2/bestTrack.ndof", float, doc="chi2/ndof")
    process.muonTable.variables.trkKink = Var("combinedQuality().trkKink", float, doc="trkKink")
    process.muonTable.variables.isStandalone = Var("isStandAloneMuon",bool,doc="muon is a standalone muon")

    process.genParticleTable.variables.vx = Var("vx", float, precision=8)
    process.genParticleTable.variables.vy = Var("vy", float, precision=8)
    process.genParticleTable.variables.vz = Var("vz", float, precision=8)

    # Skim jet variables
    process.jetTable.externalVariables = cms.PSet()
    jet_vars_rm = ['btagDeepB', 'btagCSVV2', 'qgl',  'chHEF', 'neHEF', 'chEmEF', 'neEmEF', 'muEF',
                   'chFPV0EF', 'chFPV1EF', 'chFPV2EF', 'chFPV3EF', 'hfcentralEtaStripSize', 'hfadjacentEtaStripsSize',
                   'hfsigmaEtaEta', 'hfsigmaPhiPhi']
    for v in jet_vars_rm:
        delattr(process.jetTable.variables, v)

    return process
