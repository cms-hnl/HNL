import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, Var, P3Vars
from HNL.NanoProd.common_cff import ufloat, uint, ubool

selectedDSAMuons = cms.EDFilter(
    'TrackSelector',
    src = cms.InputTag('displacedStandAloneMuons'),
    cut = cms.string('pt > 3. && abs(eta) < 2.4 && numberOfValidHits>=10')
)

diDSAMuon = cms.EDProducer(
    'DiTrackBuilder',
    src = cms.InputTag('selectedDSAMuons'),
    lep1Selection = cms.string('1'),
    lep2Selection = cms.string('1'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('userFloat("sv_ndof") > 0')
)

patDSAMuon = cms.EDProducer(
    'MuTrackBuilder',
    src1 = cms.InputTag('finalMuons'),
    src2 = cms.InputTag('selectedDSAMuons'),
    lep1Selection = cms.string('pt > 3. && isGlobalMuon && dB > 0.01'),
    lep2Selection = cms.string('1'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('charge == 0 && userFloat("sv_ndof") > 0')
)

diMuon = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('finalMuons'),
    lep1Selection = cms.string('pt > 3. && isGlobalMuon && dB > 0.01 && isLooseMuon'),
    lep2Selection = cms.string('pt > 3. && isGlobalMuon && dB > 0.01 && isLooseMuon'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('userFloat("sv_ndof") > 0')
)

dsaTable = cms.EDProducer(
    'SimpleTrackFlatTableProducer',
    src = cms.InputTag("selectedDSAMuons"),
    cut = cms.string("1"), # if we place a cut here, the indexing will be wrong
    name = cms.string("DSAMuons"),
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
    name = cms.string("DSAMuons")
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


countDiDSAMuon = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
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


diDSAMuonTables = cms.Sequence(dsaTable*diDSAMuonTable)
diDSAMuonSequence = cms.Sequence(selectedDSAMuons*diDSAMuon*dsaTable*dsaIsoTable*diDSAMuonTable)
patDSAMuonSequence = cms.Sequence(patDSAMuon*patDSAMuonTable)
diMuonSequence = cms.Sequence(diMuon*diMuonTable)


isomu24 = cms.EDFilter('TriggerResultsFilter',
    hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
    l1tResults = cms.InputTag(''),
    l1tIgnoreMaskAndPrescale = cms.bool(False),
    throw = cms.bool(True),
    triggerConditions = cms.vstring('HLT_IsoMu24_v*'),
    # mightGet = cms.optional.untracked.vstring
)

def nanoAOD_customizeDisplacedDiMuon(process):
    process.displacedDiMuonSequence = cms.Sequence(diDSAMuonSequence*patDSAMuonSequence*diMuonSequence)
    # process.nanoAOD_step.insert(1000, process.displacedDiMuonSequence)
    process.muonSequence.insert(1000, process.displacedDiMuonSequence)
    process.nanoAOD_diDSAMuon_step = cms.Path(isomu24 + process.nanoSequenceMC + countDiDSAMuon)
    process.nanoAOD_patDSAMuon_step = cms.Path(isomu24 + process.nanoSequenceMC + countPatDSAMuon)
    process.nanoAOD_diMuon_step = cms.Path(isomu24 + process.nanoSequenceMC + countDiMuon)

    process.finalMuons.cut = "pt > 3"

    # Add additional muon time variables
    process.muonTable.variables.rpcTimeInOut = Var("rpcTime().timeAtIpInOut", float, doc="RPC time in out")
    process.muonTable.variables.timeInOut = Var("time().timeAtIpInOut", float, doc="time in out")
    process.muonTable.variables.rpcTimeInOutErr = Var("rpcTime().timeAtIpInOutErr", float, doc="RPC time error in out")
    process.muonTable.variables.timeInOutErr = Var("time().timeAtIpInOutErr", float, doc="time error in out")
    process.muonTable.variables.rpcTimeNdof = Var("rpcTime().nDof", float, doc="RPC time ndof")
    process.muonTable.variables.timeNdof = Var("time().nDof", float, doc="time ndof")

    # Skim jet variables
    process.jetTable.externalVariables = cms.PSet()
    jet_vars_rm = ['btagCMVA', 'btagDeepB', 'btagCSVV2', 'qgl',  'chHEF', 'neHEF', 'chEmEF', 'neEmEF', 'muEF',
                   'chFPV0EF', 'chFPV1EF', 'chFPV2EF', 'chFPV3EF']
    for v in jet_vars_rm:
        delattr(process.jetTable.variables, v)

    return process
