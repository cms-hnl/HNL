import importlib
import os
import sys
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, Var, P3Vars
from PhysicsTools.NanoAOD.muons_cff import muonTable

def load(module_file, default_path):
  print(f'BASE PATH: {os.path.join(os.getenv("CMSSW_BASE"), "..", module_file)}')
  module_path = os.path.join(default_path, module_file)
  if not os.path.exists(module_path):
    module_path = os.path.join(os.path.dirname(__file__), module_file)
    if not os.path.exists(module_path):
      module_path = os.path.join(os.getenv("CMSSW_BASE"), 'src', module_file)
      if not os.path.exists(module_path):
        module_path = os.path.join(os.getenv("CMSSW_BASE"), '..', module_file)
        if not os.path.exists(module_path):
          raise RuntimeError(f"Cannot find path to {module_file}.")

  module_name, module_ext = os.path.splitext(module_file)
  spec = importlib.util.spec_from_file_location(module_name, module_path)
  module = importlib.util.module_from_spec(spec)
  sys.modules[module_name] = module
  spec.loader.exec_module(module)
  return module

common_cff = load('common_cff.py', os.path.join(os.getenv("CMSSW_BASE"), 'src/HNL/NanoProd/python/'))
from common_cff import ufloat, uint, ubool


def defineSelectedDSAFilter(isDSATracks, isRun2):
  # Use a placeholder for displaced pat::Muons in Run 2 
  slimmed_muons_src = cms.InputTag('slimmedMuons') if isRun2 else cms.InputTag('slimmedDisplacedMuons')
  slimmed_muons_cut = cms.string('''pt < 0''') if isRun2 else \
                      cms.string('''outerTrack.isNonnull
                              && outerTrack.pt > 3.
                              && abs(outerTrack.eta) < 2.4
                              && outerTrack.numberOfValidHits > 15
                              && outerTrack.ptError/outerTrack.pt < 1.
                              && outerTrack.chi2/outerTrack.ndof < 2.5''') 
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
          l1l2Interchangeable = cms.bool(True),
          # Use muon standalone tracks for displaced muons in Run 3
          # Only pat::Muon objects with valid outerTrack will pass the useStandalone flag
          useStandalone_l1 = cms.bool(False) if isDSATracks else cms.bool(True),
          useStandalone_l2 = cms.bool(False) if isDSATracks else cms.bool(True)
        )

def definePatDSAMuonProducer(isDSATracks):
  return cms.EDProducer(
          'MuTrackBuilder' if isDSATracks else 'DiMuonBuilder',
          src1 = cms.InputTag('finalMuons'),
          src2 = cms.InputTag('selectedDSAMuons') if isDSATracks else 
                 cms.InputTag('selectedDSAMuonsPat'),
          srcVeto = cms.InputTag('vetoMuons'),
          lep1Selection = diMuon.lepSelection,
          postVtxSelection = diDSAMuon.postVtxSelection,
          useStandalone_l1 = cms.bool(False),
          useStandalone_l2 = cms.bool(False) if isDSATracks else cms.bool(True)
        )

def defineEleDSAMuonProducer(isDSATracks):
  return cms.EDProducer(
          'EleTrackBuilder' if isDSATracks else 'EleMuBuilder',
          src1 = cms.InputTag('finalElectrons'),
          src2 = cms.InputTag('selectedDSAMuons') if isDSATracks else 
                 cms.InputTag('selectedDSAMuonsPat'),
          src2Veto = cms.InputTag('vetoMuons'),
          lep1Selection = cms.string('pt > 5. && abs(eta) < 2.5 && abs(dB("PV2D")) > 0.01'),
          postVtxSelection = diDSAMuon.postVtxSelection,
          useStandalone_l1 = cms.bool(False),
          useStandalone_l2 = cms.bool(False) if isDSATracks else cms.bool(True)
        )

def defineDSATableProducer(isDSATracks):
  prefix = '' if isDSATracks else 'outerTrack.'
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
            pt  = Var(prefix + "pt",  float, precision=-1),
            phi = Var(prefix + "phi", float, precision=12),
            eta  = Var(prefix + "eta",  float,precision=12),
            charge = Var(prefix + 'charge', int, doc='electric charge'),
            n_valid_hits = Var(prefix + 'numberOfValidHits', int, doc='valid hits'),
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
            dxy = Var(prefix + 'dxy()', float, precision=10, doc='dxy, in cm'),
            dxyErr = Var(prefix + 'dxyError()',float,doc='dxy uncertainty, in cm',precision=6),
            dz = Var(prefix + 'dz()', float, precision=10, doc='dz, in cm'),
            dzErr = Var(prefix + 'dzError()',float,doc='dz uncertainty, in cm',precision=6),
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

  # Save additional variables (only available for pat::Muon)
  this.dsaPatTable.variables.segmentComp = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14) # keep higher precision since people have cuts with 3 digits on this
  this.dsaPatTable.variables.nStations = Var("numberOfMatchedStations", "uint8", doc = "number of matched stations with default arbitration (segment & track)")
  this.dsaPatTable.variables.nTrackerLayers = Var(f"?track.isNonnull?innerTrack().hitPattern().trackerLayersWithMeasurement():0", "uint8", doc = "number of layers in the tracker")
  this.dsaPatTable.variables.highPurity = Var(f"?track.isNonnull?innerTrack().quality('highPurity'):0", bool, doc = "inner track is high purity")
  # this.dsaPatTable.variables.jetIdx = Var("?hasUserCand('jet')?userCand('jet').key():-1", "int16", doc="index of the associated jet (-1 if none)")
  this.dsaPatTable.variables.svIdx = Var("?hasUserCand('vertex')?userCand('vertex').key():-1", "int16", doc="index of matching secondary vertex")
  this.dsaPatTable.variables.tkRelIso = Var("isolationR03().sumPt/tunePMuonBestTrack().pt",float,doc="Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt",precision=6)
  this.dsaPatTable.variables.pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component")
  this.dsaPatTable.variables.pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)")
  this.dsaPatTable.variables.pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)")
  this.dsaPatTable.variables.jetRelIso = Var("?userCand('jetForLepJetVar').isNonnull()?(1./userFloat('ptRatio'))-1.:(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)",precision=8)
  # this.dsaPatTable.variables.jetPtRelv2 = Var("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0",float,doc="Relative momentum of the lepton with respect to the closest jet after subtracting the lepton",precision=8)
  this.dsaPatTable.variables.tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0", "uint8", doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)")
  this.dsaPatTable.variables.isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon")
  this.dsaPatTable.variables.isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon")
  this.dsaPatTable.variables.isStandalone = Var("isStandAloneMuon",bool,doc="muon is a standalone muon")
  this.dsaPatTable.variables.softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID")
  # this.dsaPatTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)
  this.dsaPatTable.variables.highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)")
  this.dsaPatTable.variables.pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)")
  this.dsaPatTable.variables.tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)")
  # this.dsaPatTable.variables.miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)")
  # this.dsaPatTable.variables.mvaMuID = Var("mvaIDValue()",float,doc="MVA-based ID score ",precision=6)
  # this.dsaPatTable.variables.multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)")
  # this.dsaPatTable.variables.puppiIsoId = Var("passed('PuppiIsoLoose')+passed('PuppiIsoMedium')+passed('PuppiIsoTight')", "uint8", doc="PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight)")
  this.dsaPatTable.variables.triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID")
  this.dsaPatTable.variables.inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID")
  # this.dsaPatTable.variables.jetNDauCharged = Var("?userCand('jetForLepJetVar').isNonnull()?userFloat('jetNDauChargedMVASel'):0", "uint8", doc="number of charged daughters of the closest jet")
  # time variables
  this.dsaPatTable.variables.rpcTimeInOut = Var('rpcTime().timeAtIpInOut', float, doc='RPC time in out')
  this.dsaPatTable.variables.timeInOut = Var('time().timeAtIpInOut', float, doc='time in out')
  this.dsaPatTable.variables.rpcTimeInOutErr = Var('rpcTime().timeAtIpInOutErr', float, doc='RPC time error in out')
  this.dsaPatTable.variables.timeInOutErr = Var('time().timeAtIpInOutErr', float, doc='time error in out')
  this.dsaPatTable.variables.rpcTimeNdof = Var('rpcTime().nDof', float, doc='RPC time ndof')
  this.dsaPatTable.variables.timeNdof = Var('time().nDof', float, doc='time ndof')

  # track kink
  this.dsaPatTable.variables.trkKink = Var('combinedQuality().trkKink', float, doc='trkKink')

  # For some variables, save also bestTrack information
  common_vars = ['chi2', 'ndof', 'dxy', 'dxyErr', 'dz', 'dzErr', 'pt_error', 'n_valid_hits', 'n_lost_hits']
  for var_name in common_vars:
    orig_var = getattr(this.dsaPatTable.variables, var_name)
    final_expr = ''.join(['bestTrack.'] + orig_var.expr.value().split('.')[1:])
    setattr(this.dsaPatTable.variables, f'bestTrack_{var_name}', orig_var.clone(expr=final_expr))


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
