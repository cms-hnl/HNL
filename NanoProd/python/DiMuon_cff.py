import sys
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, Var, P3Vars
from HNL.NanoProd.common_cff import ufloat, uint, ubool


def defineSelectedDSAFilter(isDSATracks, isRun2):
  # Use a placeholder for displaced pat::Muons in Run 2 
  slimmed_muons_src = cms.InputTag('slimmedMuons') if isRun2 else cms.InputTag('slimmedDisplacedMuons')
  slimmed_muons_cut = cms.string('''pt < 0''') if isRun2 else \
                      cms.string('''pt > 3. && abs(eta) < 2.4
                              && numberOfValidHits > 15
                              && bestTrack.ptError/pt < 1.
                              && bestTrack.chi2/bestTrack.ndof < 2.5''') 
  return cms.EDFilter(
          'TrackSelector' if isDSATracks else 'PATMuonRefSelector',
          src = cms.InputTag('displacedStandAloneMuons') if isDSATracks else
                slimmed_muons_src,
          cut = cms.string('''pt > 3. && abs(eta) < 2.4
                              && numberOfValidHits > 15
                              && ptError/pt < 1.
                              && chi2/ndof < 2.5''') if isDSATracks else
                slimmed_muons_cut
        )

def defineDiDSAMuonProducer(isDSATracks):
  return cms.EDProducer(
          'DiTrackBuilder' if isDSATracks else 'DiMuonBuilder',
          src = cms.InputTag('selectedDSAMuons') if isDSATracks else 
                cms.InputTag('selectedDSAMuonsPat'),
          srcVeto = cms.InputTag('vetoMuons'),
          postVtxSelection = cms.string('userFloat("sv_ndof") > 0'),
          l1l2Interchangeable = cms.bool(True)
        )

def definePatDSAMuonProducer(isDSATracks):
  return cms.EDProducer(
          'MuTrackBuilder' if isDSATracks else 'DiMuonBuilder',
          src1 = cms.InputTag('finalMuons'),
          src2 = cms.InputTag('selectedDSAMuons') if isDSATracks else 
                 cms.InputTag('selectedDSAMuonsPat'),
          srcVeto = cms.InputTag('vetoMuons'),
          lep1Selection = diMuon.lepSelection,
          postVtxSelection = diDSAMuon.postVtxSelection
        )

def defineEleDSAMuonProducer(isDSATracks):
  return cms.EDProducer(
          'EleTrackBuilder' if isDSATracks else 'EleMuBuilder',
          src1 = cms.InputTag('finalElectrons'),
          src2 = cms.InputTag('selectedDSAMuons') if isDSATracks else 
                 cms.InputTag('selectedDSAMuonsPat'),
          src2Veto = cms.InputTag('vetoMuons'),
          lep1Selection = cms.string('pt > 5. && abs(eta) < 2.5 && abs(dB("PV2D")) > 0.01'),
          postVtxSelection = diDSAMuon.postVtxSelection
        )

def defineDSATableProducer(isDSATracks):
  prefix = '' if isDSATracks else 'bestTrack().'
  return cms.EDProducer(
          'SimpleTrackFlatTableProducer' if isDSATracks else 'SimpleCandidateFlatTableProducer',
          src = cms.InputTag('selectedDSAMuons') if isDSATracks else 
                cms.InputTag('selectedDSAMuonsPat'),
          cut = cms.string('1'), # if we place a cut here, the indexing will be wrong
          name = cms.string('DSAMuon') if isDSATracks else cms.string('DSAMuonPat'),
          doc = cms.string('Displaced standalone muon variables'),
          singleton = cms.bool(False),
          extension = cms.bool(False),
          variables = cms.PSet(
            P3Vars,
            charge = Var('charge', int, doc='electric charge'),
            n_valid_hits = Var('numberOfValidHits', int, doc='valid hits'),
            n_lost_hits = Var(prefix + 'numberOfLostHits', int, doc='lost hits'),
            n_muon_stations = Var(prefix + 'hitPattern().muonStationsWithValidHits', int, doc='muon stations with valid hits'),
            n_dt_stations = Var(prefix + 'hitPattern().dtStationsWithValidHits', int, doc='DT stations with valid hits'),
            n_dt_hits = Var(prefix + 'hitPattern().numberOfValidMuonDTHits', int, doc='valid DT hits'),
            n_csc_stations = Var(prefix + 'hitPattern().cscStationsWithValidHits', int, doc='CSC stations with valid hits'),
            n_csc_hits = Var(prefix + 'hitPattern().numberOfValidMuonCSCHits', int, doc='valid CSC hits'),
            n_rpc_stations = Var(prefix + 'hitPattern().rpcStationsWithValidHits', int, doc='RPC stations with valid hits'),
            n_rpc_hits = Var(prefix + 'hitPattern().numberOfValidMuonRPCHits', int, doc='valid RPC hits'),
            chi2 = Var(prefix + 'chi2()', float, precision=10, doc='track chi2'),
            ndof = Var(prefix + 'ndof()', float,precision=10,  doc='track ndof'),
            dxy = Var(prefix + 'dxy()', float, precision=10, doc='dxy'),
            dz = Var(prefix + 'dz()', float, precision=10, doc='dz'),
            pt_error = Var(prefix + 'ptError()', float, precision=10, doc='pt error'),
            theta_error = Var(prefix + 'thetaError()', float, precision=8, doc='theta error'),
            phi_error = Var(prefix + 'phiError()', float, precision=8, doc='phi error'),
          )
        )

def defineFiltersAndProducers(isRun2):

  this = sys.modules[__name__]

  this.selectedDSAMuons = defineSelectedDSAFilter(isDSATracks=True, isRun2=isRun2)

  this.vetoMuons = cms.EDFilter(
    'PATMuonRefSelector',
    src = cms.InputTag('slimmedMuons'),
    cut = cms.string('pt > 24 && abs(eta) < 2.4 && isMediumMuon && dB < 0.02')
    # cut = cms.string('isMediumMuon && dB < 0.08')
  )

  this.diDSAMuon = defineDiDSAMuonProducer(isDSATracks=True)

  this.diSTA = cms.EDProducer(
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

  this.diMuon = cms.EDProducer(
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

  this.patDSAMuon = definePatDSAMuonProducer(isDSATracks=True)

  this.patSTA = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('finalMuons'),
    srcVeto = cms.InputTag('vetoMuons'),
    lep1Selection = diMuon.lepSelection,
    lep2Selection = diSTA.lepSelection,
    postVtxSelection = diDSAMuon.postVtxSelection
  )

  this.eleDSAMuon = defineEleDSAMuonProducer(isDSATracks=True)

  this.eleSTA = cms.EDProducer(
    'EleMuBuilder',
    src1 = cms.InputTag('finalElectrons'),
    src2 = cms.InputTag('finalMuons'),
    src2Veto = cms.InputTag('vetoMuons'),
    lep1Selection = eleDSAMuon.lep1Selection,
    lep2Selection = diSTA.lepSelection,
    postVtxSelection = diDSAMuon.postVtxSelection
  )

  this.dsaTable = defineDSATableProducer(isDSATracks=True)

  this.dsaIsoTable = cms.EDProducer(
    'TrackIsoTableProducer',
    name = cms.string('DSAMuon')
  )

  this.diDSAMuonTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('diDSAMuon'),
    cut = cms.string(''),
    name = cms.string('DiDSAMuon'),
    doc = cms.string('DiDSAMuon Variable'),
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

  # Save DSA pat:Muons 
  this.selectedDSAMuonsPat = defineSelectedDSAFilter(isDSATracks=False, isRun2=isRun2)
  this.diDSAMuonPat = defineDiDSAMuonProducer(isDSATracks=False)
  this.patDSAMuonPat = definePatDSAMuonProducer(isDSATracks=False)
  this.eleDSAMuonPat = defineEleDSAMuonProducer(isDSATracks=False)
  this.dsaPatTable = defineDSATableProducer(isDSATracks=False)

  this.dsaPatTable.variables.dxy = Var('dB("PV2D")', float, precision=10, doc='dxy (with sign) wrt first PV, in cm')
  this.dsaPatTable.variables.dz = Var('dB("PVDZ")', float, precision=10, doc='dz (with sign) wrt first PV, in cm')
  this.dsaPatTable.variables.pfIsolationR03_sumChargedHadronPt = Var('pfIsolationR03().sumChargedHadronPt()', float, doc='PF isolation dR=0.3, charged hadron component')
  this.dsaPatTable.variables.pfIsolationR03_sumChargedParticlePt = Var('pfIsolationR03().sumChargedParticlePt()', float, doc='PF isolation dR=0.3, charged particle component')
  this.dsaPatTable.variables.pfIsolationR03_sumNeutralHadronEt = Var('pfIsolationR03().sumNeutralHadronEt()', float, doc='PF isolation dR=0.3, neutral hadron component')
  this.dsaPatTable.variables.pfIsolationR03_sumPhotonEt = Var('pfIsolationR03().sumPhotonEt()', float, doc='PF isolation dR=0.3, photon component')
  this.dsaPatTable.variables.pfIsolationR03_sumPUPt = Var('pfIsolationR03().sumPUPt()', float, doc='PF isolation dR=0.3, charged PU component')
  this.dsaPatTable.variables.pfIsolationR04_sumChargedHadronPt = Var('pfIsolationR04().sumChargedHadronPt()', float, doc='PF isolation dR=0.4, charged hadron component')
  this.dsaPatTable.variables.pfIsolationR04_sumChargedParticlePt = Var('pfIsolationR04().sumChargedParticlePt()', float, doc='PF isolation dR=0.4, charged particle component')
  this.dsaPatTable.variables.pfIsolationR04_sumNeutralHadronEt = Var('pfIsolationR04().sumNeutralHadronEt()', float, doc='PF isolation dR=0.4, neutral hadron component')
  this.dsaPatTable.variables.pfIsolationR04_sumPhotonEt = Var('pfIsolationR04().sumPhotonEt()', float, doc='PF isolation dR=0.4, photon component')
  this.dsaPatTable.variables.pfIsolationR04_sumPUPt = Var('pfIsolationR04().sumPUPt()', float, doc='PF isolation dR=0.4, charged PU component')
  this.dsaPatTable.variables.rpcTimeInOut = Var('rpcTime().timeAtIpInOut', float, doc='RPC time in out')
  this.dsaPatTable.variables.timeInOut = Var('time().timeAtIpInOut', float, doc='time in out')
  this.dsaPatTable.variables.rpcTimeInOutErr = Var('rpcTime().timeAtIpInOutErr', float, doc='RPC time error in out')
  this.dsaPatTable.variables.timeInOutErr = Var('time().timeAtIpInOutErr', float, doc='time error in out')
  this.dsaPatTable.variables.rpcTimeNdof = Var('rpcTime().nDof', float, doc='RPC time ndof')
  this.dsaPatTable.variables.timeNdof = Var('time().nDof', float, doc='time ndof')


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
    'drop abs(pdgId) == 2212 && abs(pz) > 1000', #drop LHC protons accidentally added by previous keeps
  ]

  for coord in [ 'x', 'y', 'z' ]:
    setattr(process.genParticleTable.variables, 'v' + coord,
            Var(f'vertex().{coord}', float, precision=10,
                doc=f'{coord} coordinate of the gen particle production vertex'))

  return process


def nanoAOD_customizeDisplacedDiMuon(process, isRun2):

  defineFiltersAndProducers(isRun2)
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

  process.MessageLogger.cerr.FwkReport.reportEvery = 100

  process.displacedDiMuonSequence = cms.Sequence()
  n_mod = len(_all_modules)
  for name in _all_modules:
    x = getattr(this, name)
    setattr(process, name, x)
    process.displacedDiMuonSequence.insert(n_mod, x)

  process.nanoSequenceCommon.insert(1000, process.displacedDiMuonSequence)

  process.finalMuons.cut = 'pt > 3'

  # Add additional muon time variables
  process.muonTable.variables.rpcTimeInOut = Var('rpcTime().timeAtIpInOut', float, doc='RPC time in out')
  process.muonTable.variables.timeInOut = Var('time().timeAtIpInOut', float, doc='time in out')
  process.muonTable.variables.rpcTimeInOutErr = Var('rpcTime().timeAtIpInOutErr', float, doc='RPC time error in out')
  process.muonTable.variables.timeInOutErr = Var('time().timeAtIpInOutErr', float, doc='time error in out')
  process.muonTable.variables.rpcTimeNdof = Var('rpcTime().nDof', float, doc='RPC time ndof')
  process.muonTable.variables.timeNdof = Var('time().nDof', float, doc='time ndof')
  process.muonTable.variables.nValidHits = Var('numberOfValidHits', float, doc='n valid hits')
  process.muonTable.variables.chi2ndof = Var('bestTrack.chi2/bestTrack.ndof', float, doc='chi2/ndof')
  process.muonTable.variables.trkKink = Var('combinedQuality().trkKink', float, doc='trkKink')
  process.muonTable.variables.isStandalone = Var('isStandAloneMuon',bool,doc='muon is a standalone muon')

  process.electronTable.variables.dEtaSeed = Var(
    'deltaEtaSuperClusterTrackAtVtx - superCluster.eta() + superCluster.seed().eta()', float, doc='', precision=10)
  process.electronTable.variables.dPhiIn = Var('deltaPhiSuperClusterTrackAtVtx', float, doc='', precision=10)

  process = customiseGenParticles(process)

  return process


def nanoAOD_customizeDisplacedDiMuon_Run2(process):
  return nanoAOD_customizeDisplacedDiMuon(process, True)


def nanoAOD_customizeDisplacedDiMuon_Run3(process):
  return nanoAOD_customizeDisplacedDiMuon(process, False)
