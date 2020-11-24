import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, Var, P3Vars
from HNL.NanoProd.common_cff import ufloat, uint, ubool

diDSAMuon = cms.EDProducer(
    'DiTrackBuilder',
    src = cms.InputTag('displacedStandAloneMuons'),
    lep1Selection = cms.string('pt > 3. && abs(eta) < 2.4 && numberOfValidHits>=10'),
    lep2Selection = cms.string('pt > 3. && abs(eta) < 2.4 && numberOfValidHits>=10'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('1')
)

patDSAMuon = cms.EDProducer(
    'MuTrackBuilder',
    src1 = cms.InputTag('slimmedMuons'),
    src2 = cms.InputTag('displacedStandAloneMuons'),
    lep1Selection = cms.string('1'),
    lep2Selection = cms.string('pt > 3. && abs(eta) < 2.4 && numberOfValidHits>=10'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('1')
)

diMuon = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('slimmedMuons'),
    lep1Selection = cms.string('1'),
    lep2Selection = cms.string('1'),
    preVtxSelection = cms.string('1'),
    postVtxSelection = cms.string('1')
)

dsaTable = cms.EDProducer(
    'SimpleTrackFlatTableProducer',
    src = cms.InputTag("displacedStandAloneMuons"),
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

CountDisplacedDiMuon = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("diDSAMuon")
)    


diDSAMuonTables = cms.Sequence(dsaTable*diDSAMuonTable)
diDSAMuonSequence = cms.Sequence(diDSAMuon*dsaTable*diDSAMuonTable)

def nanoAOD_customizeDisplacedDiMuon(process):
    process.diDSAMuonSequence = cms.Sequence(diDSAMuonSequence)
    process.nanoAOD_step.insert(0, process.diDSAMuonSequence)
    process.nanoAOD_displacedDiMuon_step = cms.Path(process.nanoSequenceMC + process.diDSAMuonSequence + CountDisplacedDiMuon)
    return process
