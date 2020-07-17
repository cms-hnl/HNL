import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars
from HNL.NanoProd.common_cff import ufloat, uint, ubool

displacedTracks = cms.EDProducer(
    'DiTrackBuilder',
    src = cms.InputTag('displacedStandAloneMuons'),
    lep1Selection = cms.string('pt > 3. && abs(eta) < 2.4'),
    lep2Selection = cms.string('pt > 3. && abs(eta) < 2.4'),
    preVtxSelection = cms.string(
        '1'
    ),
    postVtxSelection = cms.string('1'),
)
displacedDiMuonTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("displacedTracks"),
    cut = cms.string(""),
    name = cms.string("DisplacedDiMuon"),
    doc = cms.string("DisplacedDiMuon Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),

        # fit and vtx info
        chi2 = ufloat('sv_chi2'),
        ndof = ufloat('sv_ndof'),
        svprob = ufloat('sv_prob'),
        # l_xy = ufloat('l_xy'),
        # l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll
        # mll_raw = Var('userCand("dilepton").mass()', float),
        # mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        # mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work
        # mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        # cos2D = ufloat('cos_theta_2D'),
        # fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        # fit_pt = ufloat('fitted_pt'),
        # fit_eta = ufloat('fitted_eta'),
        # fit_phi = ufloat('fitted_phi'),
        # fit_l1_pt = ufloat('fitted_l1_pt'),
        # fit_l1_eta = ufloat('fitted_l1_eta'),
        # fit_l1_phi = ufloat('fitted_l1_phi'),
        # fit_l2_pt = ufloat('fitted_l2_pt'),
        # fit_l2_eta = ufloat('fitted_l2_eta'),
        # fit_l2_phi = ufloat('fitted_l2_phi'),
        # fit_k_pt = ufloat('fitted_k_pt'),
        # fit_k_eta = ufloat('fitted_k_eta'),
        # fit_k_phi = ufloat('fitted_k_phi'),
    )
)

CountDisplacedDiMuon = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("displacedTracks")
)    


displaceDiMuonTables = cms.Sequence(displacedDiMuonTable)
displacedDiMuonSequence = cms.Sequence(displacedTracks*displacedDiMuonTable)

def nanoAOD_customizeDisplacedDiMuon(process):
    process.displacedDiMuonSequence = cms.Sequence(displacedDiMuonSequence)
    process.nanoAOD_step.insert(0, process.displacedDiMuonSequence)
    # process.nanoAOD_displacedDiMuon_step = cms.Path(process.nanoSequence + process.displacedDiMuonSequence + CountDisplacedDiMuon )
    return process
