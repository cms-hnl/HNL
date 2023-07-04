import sys
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, Var, P3Vars
from HNL.NanoProd.common_cff import ufloat, uint, ubool


selectedDSAMuons = cms.EDFilter(
  'PATMuonRefSelector',
  src = cms.InputTag('slimmedDisplacedMuons'),
  cut = cms.string('''pt > 3. && abs(eta) < 2.4
                      && numberOfValidHits > 15
                      && bestTrack.ptError/pt < 1.
                      && bestTrack.chi2/bestTrack.ndof < 2.5''')
)

vetoMuons = cms.EDFilter(
  'PATMuonRefSelector',
  src = cms.InputTag('slimmedMuons'),
  cut = cms.string('pt > 24 && abs(eta) < 2.4 && isMediumMuon && dB < 0.02')
  # cut = cms.string('isMediumMuon && dB < 0.08')
)

diDSAMuon = cms.EDProducer(
  'DiMuonBuilder',
  src = cms.InputTag('selectedDSAMuons'),
  srcVeto = cms.InputTag('vetoMuons'),
  postVtxSelection = cms.string('userFloat("sv_ndof") > 0'),
  l1l2Interchangeable = cms.bool(True)
)

diSTA = cms.EDProducer(
  'DiMuonBuilder',
  src = cms.InputTag('finalMuons'),
  srcVeto = cms.InputTag('vetoMuons'),
  lepSelection = cms.string('''pt > 3. && abs(eta) < 2.4
                               && isStandAloneMuon
                               && numberOfValidHits > 12
                               && bestTrack.ptError/pt < 1.
                               && bestTrack.chi2/bestTrack.ndof < 2.5'''),
  postVtxSelection = diDSAMuon.postVtxSelection,
  l1l2Interchangeable = cms.bool(True)
)

diMuon = cms.EDProducer(
  'DiMuonBuilder',
  src = cms.InputTag('finalMuons'),
  srcVeto = cms.InputTag('vetoMuons'),
  lepSelection = cms.string('''pt > 3. && abs(eta) < 2.4
                               && isGlobalMuon
                               && isLooseMuon
                               && dB > 0.01
                               && segmentCompatibility > 0.451
                               && combinedQuality().trkKink < 20'''),
  postVtxSelection = diDSAMuon.postVtxSelection,
  l1l2Interchangeable = cms.bool(True)
)

patDSAMuon = cms.EDProducer(
  'DiMuonBuilder',
  src1 = cms.InputTag('finalMuons'),
  src2 = cms.InputTag('selectedDSAMuons'),
  srcVeto = cms.InputTag('vetoMuons'),
  lep1Selection = diMuon.lepSelection,
  postVtxSelection = diDSAMuon.postVtxSelection
)

patSTA = cms.EDProducer(
  'DiMuonBuilder',
  src = cms.InputTag('finalMuons'),
  srcVeto = cms.InputTag('vetoMuons'),
  lep1Selection = diMuon.lepSelection,
  lep2Selection = diSTA.lepSelection,
  postVtxSelection = diDSAMuon.postVtxSelection
)

eleDSAMuon = cms.EDProducer(
  'EleMuBuilder',
  src1 = cms.InputTag('finalElectrons'),
  src2 = cms.InputTag('selectedDSAMuons'),
  src2Veto = cms.InputTag('vetoMuons'),
  lep1Selection = cms.string("pt > 5. && abs(eta) < 2.5 && abs(dB('PV2D')) > 0.01"),
  postVtxSelection = diDSAMuon.postVtxSelection
)

eleSTA = cms.EDProducer(
  'EleMuBuilder',
  src1 = cms.InputTag('finalElectrons'),
  src2 = cms.InputTag('finalMuons'),
  src2Veto = cms.InputTag('vetoMuons'),
  lep1Selection = eleDSAMuon.lep1Selection,
  lep2Selection = diSTA.lepSelection,
  postVtxSelection = diDSAMuon.postVtxSelection
)

### TO DO: fix expression for all variables and add isolation for DSA muons ###
dsaTable = cms.EDProducer(
  'SimpleCandidateFlatTableProducer',
  src = cms.InputTag("selectedDSAMuons"),
  cut = cms.string("1"), # if we place a cut here, the indexing will be wrong
  name = cms.string("DSAMuon"),
  doc = cms.string("Displaced standalone muon variables"),
  singleton = cms.bool(False),
  extension = cms.bool(False),
  variables = cms.PSet(
    P3Vars,
    charge = Var("charge", int, doc="electric charge"),
    n_valid_hits = Var('numberOfValidHits', int, doc='valid hits'),
    n_lost_hits = Var('bestTrack().numberOfLostHits()', int, doc='lost hits'),
    n_muon_stations = Var('bestTrack().hitPattern().muonStationsWithValidHits', int, doc='muon stations with valid hits'),
    n_dt_stations = Var('bestTrack().hitPattern().dtStationsWithValidHits', int, doc='DT stations with valid hits'),
    n_dt_hits = Var('bestTrack().hitPattern().numberOfValidMuonDTHits', int, doc='valid DT hits'),
    n_csc_stations = Var('bestTrack().hitPattern().cscStationsWithValidHits', int, doc='CSC stations with valid hits'),
    n_csc_hits = Var('bestTrack().hitPattern().numberOfValidMuonCSCHits', int, doc='valid CSC hits'),
    n_rpc_stations = Var('bestTrack().hitPattern().rpcStationsWithValidHits', int, doc='RPC stations with valid hits'),
    n_rpc_hits = Var('bestTrack().hitPattern().numberOfValidMuonRPCHits', int, doc='valid RPC hits'),
    chi2 = Var('bestTrack().chi2()', float, precision=10, doc='track chi2'),
    ndof = Var('bestTrack().ndof()', float,precision=10,  doc='track ndof'),
    dxy = Var("dB('PV2D')", float, precision=10, doc='dxy (with sign) wrt first PV, in cm'),
    dz = Var("dB('PVDZ')", float, precision=10, doc='dz (with sign) wrt first PV, in cm'),
    pt_error = Var('bestTrack().ptError()', float, precision=10, doc='pt error'),
    theta_error = Var('bestTrack().thetaError()', float, precision=8, doc='theta error'),
    phi_error = Var('bestTrack().phiError()', float, precision=8, doc='phi error'),
    pfIsolationR03_sumChargedHadronPt = Var('pfIsolationR03().sumChargedHadronPt()', float, doc='PF isolation dR=0.3, charged hadron component'),
    pfIsolationR03_sumChargedParticlePt = Var('pfIsolationR03().sumChargedParticlePt()', float, doc='PF isolation dR=0.3, charged particle component'),
    pfIsolationR03_sumNeutralHadronEt = Var('pfIsolationR03().sumNeutralHadronEt()', float, doc='PF isolation dR=0.3, neutral hadron component'),
    pfIsolationR03_sumPhotonEt = Var('pfIsolationR03().sumPhotonEt()', float, doc='PF isolation dR=0.3, photon component'),
    pfIsolationR03_sumPUPt = Var('pfIsolationR03().sumPUPt()', float, doc='PF isolation dR=0.3, charged PU component'),
    pfIsolationR04_sumChargedHadronPt = Var('pfIsolationR04().sumChargedHadronPt()', float, doc='PF isolation dR=0.4, charged hadron component'),
    pfIsolationR04_sumChargedParticlePt = Var('pfIsolationR04().sumChargedParticlePt()', float, doc='PF isolation dR=0.4, charged particle component'),
    pfIsolationR04_sumNeutralHadronEt = Var('pfIsolationR04().sumNeutralHadronEt()', float, doc='PF isolation dR=0.4, neutral hadron component'),
    pfIsolationR04_sumPhotonEt = Var('pfIsolationR04().sumPhotonEt()', float, doc='PF isolation dR=0.4, photon component'),
    pfIsolationR04_sumPUPt = Var('pfIsolationR04().sumPUPt()', float, doc='PF isolation dR=0.4, charged PU component'),
    rpcTimeInOut = Var("rpcTime().timeAtIpInOut", float, doc="RPC time in out"),
    timeInOut = Var("time().timeAtIpInOut", float, doc="time in out"),
    rpcTimeInOutErr = Var("rpcTime().timeAtIpInOutErr", float, doc="RPC time error in out"),
    timeInOutErr = Var("time().timeAtIpInOutErr", float, doc="time error in out"),
    rpcTimeNdof = Var("rpcTime().nDof", float, doc="RPC time ndof"),
    timeNdof = Var("time().nDof", float, doc="time ndof"),
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

this = sys.modules[__name__]
_all_filters = []
_all_producers = []
_all_tables = []
for x_name in dir(this):
  x = getattr(this, x_name)
  x_type = type(x)
  if x_type == cms.EDFilter:
    _all_filters.append(x_name)
  if x_type == cms.EDProducer:
    _all_producers.append(x_name)
    if x.type_().endswith('Builder'):
      if x_name != 'diDSAMuon':
        setattr(this, x_name + 'Table', diDSAMuonTable.clone(
          src = x_name,
          name = x_name[0].capitalize() + x_name[1:],
          doc = x_name + ' variables'
        ))
      _all_tables.append(x_name + 'Table')
_all_modules = _all_filters + _all_producers + _all_tables

def customiseGenParticles(process):
  def pdgOR(pdgs):
    abs_pdgs = [ f'abs(pdgId) == {pdg}' for pdg in pdgs ]
    return '( ' + ' || '.join(abs_pdgs) + ' )'

  leptons = pdgOR([ 11, 13, 15 ])
  important_particles = pdgOR([ 6, 23, 24, 25, 35, 39, 9990012, 9900012 ])
  process.finalGenParticles.select = [
    'drop *',
    'keep++ statusFlags().isLastCopy() && ' + leptons,
    '+keep statusFlags().isFirstCopy() && ' + leptons,
    'keep+ statusFlags().isLastCopy() && ' + important_particles,
    '+keep statusFlags().isFirstCopy() && ' + important_particles,
    "drop abs(pdgId) == 2212 && abs(pz) > 1000", #drop LHC protons accidentally added by previous keeps
  ]

  for coord in [ 'x', 'y', 'z' ]:
    setattr(process.genParticleTable.variables, 'v' + coord,
            Var(f'vertex().{coord}', float, precision=10,
                doc=f'{coord} coordinate of the gen particle production vertex'))

  return process

def nanoAOD_customizeDisplacedDiMuon(process):

  process.MessageLogger.cerr.FwkReport.reportEvery = 100

  process.displacedDiMuonSequence = cms.Sequence()
  n_mod = len(_all_modules)
  for name in _all_modules:
    x = getattr(this, name)
    setattr(process, name, x)
    process.displacedDiMuonSequence.insert(n_mod, x)

  process.nanoSequenceCommon.insert(1000, process.displacedDiMuonSequence)

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

  process.electronTable.variables.dEtaSeed = Var(
    "deltaEtaSuperClusterTrackAtVtx - superCluster.eta() + superCluster.seed().eta()", float, doc="", precision=10)
  process.electronTable.variables.dPhiIn = Var("deltaPhiSuperClusterTrackAtVtx", float, doc="", precision=10)

  process = customiseGenParticles(process)

  # process.tauTask.remove(process.rerunMvaIsolationTaskForNano)
  # del process.slimmedTausUpdated.tauIDSources.byDeepTau2018v2p5VSeraw
  # del process.slimmedTausUpdated.tauIDSources.byDeepTau2018v2p5VSmuraw
  # del process.slimmedTausUpdated.tauIDSources.byDeepTau2018v2p5VSjetraw
  # del process.tauTable.variables.idDeepTau2018v2p5VSe
  # del process.tauTable.variables.idDeepTau2018v2p5VSmu
  # del process.tauTable.variables.idDeepTau2018v2p5VSjet
  # del process.tauTable.variables.rawDeepTau2018v2p5VSe
  # del process.tauTable.variables.rawDeepTau2018v2p5VSmu
  # del process.tauTable.variables.rawDeepTau2018v2p5VSjet
  # process.finalTaus.cut = 'pt < 0'

  # for task_name in [ 'boostedTauTablesTask', 'boostedTauTask', 'jetAK8LepTask', 'jetAK8TablesTask', 'jetAK8Task' ]:
  #   task = getattr(process, task_name)
  #   process.nanoTableTaskCommon.remove(task)

  # process.slimmedElectronsUpdated.src = "slimmedElectronsTo106X"
  # process.finalLowPtElectrons.cut = 'pt < 0'
  # process.finalLowPtElectrons.src = 'slimmedElectrons'
  # del process.jetTable.variables.hfadjacentEtaStripsSize
  # del process.jetTable.variables.hfcentralEtaStripSize
  # del process.jetTable.variables.hfsigmaEtaEta
  # del process.jetTable.variables.hfsigmaPhiPhi

  # process.jetMCTask.remove(process.jetMCTaskak8)
  # process.nanoTableTaskFS.remove(process.boostedTauMCTask)
  # #process.nanoTableTaskFS.remove(process.lowPtElectronMCTask)
  return process
