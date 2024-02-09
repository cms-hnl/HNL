#include <vector>
#include <memory>
#include <map>
#include <string>
#include <limits>
#include <algorithm>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "KinVtxFitter.h"

namespace {
  template<typename Lepton>
  const reco::Track& getTrack(const Lepton& lep) { return *lep.bestTrack(); }
  template<>
  const reco::Track& getTrack<reco::Track>(const reco::Track& lep) { return lep; }

  template<typename Lepton>
  const reco::Track& getOuterTrack(const Lepton& lep) { return *lep.bestTrack(); }
  template<>
  const reco::Track& getOuterTrack<reco::Track>(const reco::Track& lep) { return lep; }
  template<>
  const reco::Track& getOuterTrack<pat::Muon>(const pat::Muon& mu) { return *mu.outerTrack().get(); }

  template<typename Lepton>
  double getMass();

  template<>
  double getMass<pat::Electron>() { return 0.00051099895000; }

  template<>
  double getMass<pat::Muon>() { return 0.1056583755; }

  template<>
  double getMass<reco::Track>() { return getMass<pat::Muon>(); /* assuming muon */ }

  template<typename Lepton>
  double getMassUnc();

  template<>
  double getMassUnc<pat::Electron>() { return 1.5e-13; }

  template<>
  double getMassUnc<pat::Muon>() { return 2.3e-9; }

  template<>
  double getMassUnc<reco::Track>() { return getMassUnc<pat::Muon>(); /* assuming muon */ }


  template<typename Lepton>
  math::XYZTLorentzVector getP4(const Lepton& lep) { return lep.p4(); }
  template<>
  math::XYZTLorentzVector getP4<reco::Track>(const reco::Track& lep)
  {
    const math::PtEtaPhiMLorentzVector p4(lep.pt(), lep.eta(), lep.phi(), getMass<reco::Track>());
    return math::XYZTLorentzVector(p4);
  }

  template<typename Lepton>
  bool dR_match(const Lepton& lep, const edm::View<reco::Candidate>& veto_leptons, double deltaR_thr)
  {
    for(const auto& cand : veto_leptons) {
      if(reco::deltaR(lep, cand) < deltaR_thr) {
        return true;
      }
    }
    return false;
  }

  template<typename Lepton>
  bool outerTrackAvailable(const Lepton& lep) { return false; }

  template<>
  bool outerTrackAvailable<pat::Muon>(const pat::Muon& mu) { return mu.outerTrack().isNonnull(); }

  template<typename Lepton1, typename Lepton2>
  struct DiLeptonBuilderHelper {
    static bool isSameObject(const Lepton1& lepton1, const Lepton2& lepton2) { return false; }
  };

  template<typename Lepton>
  struct DiLeptonBuilderHelper<Lepton, Lepton> {
    static bool isSameObject(const Lepton& lepton1, const Lepton& lepton2) { return &lepton1 == &lepton2; }
  };
}

template<typename Lepton1, typename Lepton2>
class DiLeptonBuilder : public edm::global::EDProducer<> {
public:
  using Helper = DiLeptonBuilderHelper<Lepton1, Lepton2>;
  using Lepton1Collection = edm::View<Lepton1>;
  using Lepton2Collection = edm::View<Lepton2>;
  static constexpr bool l1l2HaveSameType = std::is_same<Lepton1, Lepton2>::value;

  explicit DiLeptonBuilder(const edm::ParameterSet &cfg) :
    ttkToken_(esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"}))
  {
    if(cfg.exists("lepSelection")) {
      if(cfg.exists("lep1Selection") || cfg.exists("lep2Selection"))
        throw std::runtime_error("Inconsistent config. lepSelection and lepNSelection are specified at the same time.");
      const std::string sel = cfg.getParameter<std::string>("lepSelection");
      l1_selection_ = std::make_unique<StringCutObjectSelector<Lepton1>>(sel);
      l2_selection_ = std::make_unique<StringCutObjectSelector<Lepton2>>(sel);
    }
    if(cfg.exists("lep1Selection")) {
      const std::string sel = cfg.getParameter<std::string>("lep1Selection");
      l1_selection_ = std::make_unique<StringCutObjectSelector<Lepton1>>(sel);
    }
    if(cfg.exists("lep2Selection")) {
      const std::string sel = cfg.getParameter<std::string>("lep2Selection");
      l2_selection_ = std::make_unique<StringCutObjectSelector<Lepton2>>(sel);
    }

    if(cfg.exists("preVtxSelection")) {
      const std::string sel = cfg.getParameter<std::string>("preVtxSelection");
      pre_vtx_selection_ = std::make_unique<StringCutObjectSelector<pat::CompositeCandidate>>(sel);
    }
    if(cfg.exists("postVtxSelection")) {
      const std::string sel = cfg.getParameter<std::string>("postVtxSelection");
      post_vtx_selection_ = std::make_unique<StringCutObjectSelector<pat::CompositeCandidate>>(sel);
    }

    if(cfg.exists("src")) {
      if(cfg.exists("src1") || cfg.exists("src2"))
        throw std::runtime_error("Inconsistent config. src and srcN are specified at the same time.");
      if(!l1l2HaveSameType)
        throw std::runtime_error("Inconsistent config. Only one src is specified, but Lepton1 type != Lepton2 type.");
      src1_ = std::make_unique<edm::EDGetTokenT<Lepton1Collection>>(
              consumes<Lepton1Collection>(cfg.getParameter<edm::InputTag>("src")));
      src2_ = std::make_unique<edm::EDGetTokenT<Lepton2Collection>>(
              consumes<Lepton2Collection>(cfg.getParameter<edm::InputTag>("src")));
      commonSrc_ = true;
    } else {
      src1_ = std::make_unique<edm::EDGetTokenT<Lepton1Collection>>(
              consumes<Lepton1Collection>(cfg.getParameter<edm::InputTag>("src1")));
      src2_ = std::make_unique<edm::EDGetTokenT<Lepton2Collection>>(
              consumes<Lepton2Collection>(cfg.getParameter<edm::InputTag>("src2")));
      commonSrc_ = false;
    }

    if(cfg.exists("srcVeto")) {
      if(cfg.exists("src1Veto") || cfg.exists("src2Veto"))
        throw std::runtime_error("Inconsistent config. srcVeto and srcNVeto are specified at the same time.");

      src1Veto_ = std::make_unique<edm::EDGetTokenT<edm::View<reco::Candidate>>>(
                  consumes<edm::View<reco::Candidate>>( cfg.getParameter<edm::InputTag>("srcVeto")));
      src2Veto_ = std::make_unique<edm::EDGetTokenT<edm::View<reco::Candidate>>>(
                  consumes<edm::View<reco::Candidate>>( cfg.getParameter<edm::InputTag>("srcVeto")));
    } else {
      if(cfg.exists("src1Veto")) {
        src1Veto_ = std::make_unique<edm::EDGetTokenT<edm::View<reco::Candidate>>>(
                    consumes<edm::View<reco::Candidate>>( cfg.getParameter<edm::InputTag>("src1Veto")));
      }
      if(cfg.exists("src2Veto")) {
        src2Veto_ = std::make_unique<edm::EDGetTokenT<edm::View<reco::Candidate>>>(
                    consumes<edm::View<reco::Candidate>>( cfg.getParameter<edm::InputTag>("src2Veto")));
      }
    }

    if(cfg.exists("deltaR_thr")) {
      deltaR_thr_ = cfg.getParameter<double>("deltaR_thr");
    }

    if(cfg.exists("l1l2Interchangeable")) {
      l1l2_interchangeable_ = cfg.getParameter<bool>("l1l2Interchangeable");
      if(l1l2_interchangeable_ && !l1l2HaveSameType)
        throw std::runtime_error("Inconsistent config. l1l2Interchangeable = true, but Lepton1 type != Lepton2 type.");
    }

    if(cfg.exists("verbose")) {
      verbose_ = cfg.getParameter<int>("verbose");
    }

    if(cfg.exists("useStandalone_l1")) {
      useStandalone_l1_ = cfg.getParameter<bool>("useStandalone_l1");
    }
    if(cfg.exists("useStandalone_l2")) {
      useStandalone_l2_ = cfg.getParameter<bool>("useStandalone_l2");
    }
    if(l1l2_interchangeable_ && (useStandalone_l1_ != useStandalone_l2_))
      throw std::runtime_error("Inconsistent config. l1l2Interchangeable = true, but useStandalone_l1 != useStandalone_l2.");

    produces<pat::CompositeCandidateCollection>();
  }

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:
  std::unique_ptr<StringCutObjectSelector<Lepton1>> l1_selection_; // cut on leading lepton
  std::unique_ptr<StringCutObjectSelector<Lepton2>> l2_selection_; // cut on sub-leading lepton
  std::unique_ptr<StringCutObjectSelector<pat::CompositeCandidate>> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  std::unique_ptr<StringCutObjectSelector<pat::CompositeCandidate>> post_vtx_selection_; // cut on the di-lepton after the SV fit
  std::unique_ptr<edm::EDGetTokenT<Lepton1Collection>> src1_;
  std::unique_ptr<edm::EDGetTokenT<Lepton2Collection>> src2_;
  std::unique_ptr<edm::EDGetTokenT<edm::View<reco::Candidate>>> src1Veto_;
  std::unique_ptr<edm::EDGetTokenT<edm::View<reco::Candidate>>> src2Veto_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  bool l1l2_interchangeable_{false};
  bool commonSrc_{false};
  double deltaR_thr_{0.3};
  int verbose_{0};
  bool useStandalone_l1_{false};
  bool useStandalone_l2_{false};
};


template<typename Lepton1, typename Lepton2>
void DiLeptonBuilder<Lepton1, Lepton2>::produce(edm::StreamID, edm::Event& evt, edm::EventSetup const& setup) const {

  const TransientTrackBuilder& tt_builder = setup.getData(ttkToken_);

  edm::Handle<Lepton1Collection> leptons1;
  evt.getByToken(*src1_, leptons1);
  edm::Handle<Lepton2Collection> leptons2;
  evt.getByToken(*src2_, leptons2);
  edm::Handle<edm::View<reco::Candidate>> vetoLeptons1;
  if(src1Veto_)
    evt.getByToken(*src1Veto_, vetoLeptons1);
  edm::Handle<edm::View<reco::Candidate>> vetoLeptons2;
  if(src2Veto_)
    evt.getByToken(*src2Veto_, vetoLeptons2);

  auto dilepton_pairs = std::make_unique<pat::CompositeCandidateCollection>();

  for(size_t l1_idx = 0; l1_idx < leptons1->size(); ++l1_idx) {
    const Lepton1& l1 = leptons1->at(l1_idx);
    if(l1_selection_ && !(*l1_selection_)(l1)) continue;

    // Get track associated to l1
    reco::Track track_l1;
    if (useStandalone_l1_) {
      // Accept only pat::Muon objects with an outer track
      if (!(outerTrackAvailable)(l1)) continue; 
      track_l1 = getOuterTrack(l1);
    } else {
    track_l1 = getTrack(l1);
    }
    // Apply veto on l1 using deltaR matching with chosen l1 track 
    if(src1Veto_ && dR_match(track_l1, *vetoLeptons1, deltaR_thr_)) continue;

    const size_t l2_start = l1l2_interchangeable_ && commonSrc_ ? l1_idx + 1 : 0;
    for(size_t l2_idx = l2_start; l2_idx < leptons2->size(); ++l2_idx) {
      const Lepton2& l2 = leptons2->at(l2_idx);
      if(Helper::isSameObject(l1, l2)) continue;
      if(l2_selection_ && !(*l2_selection_)(l2)) continue;

      // Get track associated to l2
      reco::Track track_l2;
      if (useStandalone_l2_) {
        if (!(outerTrackAvailable)(l2)) continue;
        track_l2 = getOuterTrack(l2);
      } else {
        track_l2 = getTrack(l2);
      }
      // Apply veto on l2 using deltaR matching with chosen l2 track 
      if(src2Veto_ && dR_match(track_l2, *vetoLeptons2, deltaR_thr_)) continue;

      pat::CompositeCandidate lepton_pair;
      /* Can't always use track here because muon mass is assumed in P4 computation of reco::Track */
      if (useStandalone_l1_ && useStandalone_l2_) {lepton_pair.setP4(getP4(track_l1) + getP4(track_l2));}
      else if (useStandalone_l1_) {lepton_pair.setP4(getP4(track_l1) + getP4(l2));}
      else if (useStandalone_l2_) {lepton_pair.setP4(getP4(l1) + getP4(track_l2));}
      else {lepton_pair.setP4(getP4(l1) + getP4(l2));}

      lepton_pair.setCharge(track_l1.charge() + track_l2.charge());
      lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(track_l1, track_l2));

      // Put the lepton passing the corresponding selection
      if(!l1l2_interchangeable_ || track_l1.pt() >= track_l2.pt()) {
        lepton_pair.addUserInt("l1_idx", l1_idx);
        lepton_pair.addUserInt("l2_idx", l2_idx);
      } else {
        lepton_pair.addUserInt("l1_idx", l2_idx);
        lepton_pair.addUserInt("l2_idx", l1_idx);
      }

      // before making the SV, cut on the info we have
      if(pre_vtx_selection_ && !(*pre_vtx_selection_)(lepton_pair) ) continue;

      reco::TransientTrack tt_l1(tt_builder.build(track_l1));
      reco::TransientTrack tt_l2(tt_builder.build(track_l2));

      try {
        KinVtxFitter fitter(
          {tt_l1, tt_l2},
          {getMass<Lepton1>(), getMass<Lepton2>()},
          {getMassUnc<Lepton1>(), getMassUnc<Lepton2>()} //some small sigma for the particle mass
        );
        reco::Candidate::Point vtx;
        if (fitter.success()) {
          const auto& fitted_vtx = fitter.fitted_vtx();
          vtx = reco::Candidate::Point(fitted_vtx.x(), fitted_vtx.y(), fitted_vtx.z());
        }
        lepton_pair.setVertex(vtx);

        lepton_pair.addUserFloat("sv_chi2", fitter.chi2());
        lepton_pair.addUserFloat("sv_ndof", fitter.dof()); // float??
        lepton_pair.addUserFloat("sv_prob", fitter.prob());
        lepton_pair.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1.);
        lepton_pair.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1.);
        lepton_pair.addUserFloat("fitted_pt", fitter.success() ? sqrt(pow(fitter.fitted_candidate().globalMomentum().x(), 2) +  pow(fitter.fitted_candidate().globalMomentum().y(), 2)) : -1.);
        lepton_pair.addUserFloat("vtx_x", lepton_pair.vx());
        lepton_pair.addUserFloat("vtx_y", lepton_pair.vy());
        lepton_pair.addUserFloat("vtx_z", lepton_pair.vz());
        lepton_pair.addUserFloat("vtx_ex", fitter.success() ? sqrt(fitter.fitted_vtx_uncertainty().cxx()) : -1.);
        lepton_pair.addUserFloat("vtx_ey", fitter.success() ? sqrt(fitter.fitted_vtx_uncertainty().cyy()) : -1.);
        lepton_pair.addUserFloat("vtx_ez", fitter.success() ? sqrt(fitter.fitted_vtx_uncertainty().czz()) : -1.);
      } catch (const std::exception& e) {
        if(verbose_ > 0) {
          std::cerr << e.what() << std::endl;
          std::cerr << "l1 pt, eta, phi, dxy, dz " << track_l1.pt() << ", " << track_l1.eta() << ", " << track_l1.phi()
                    << ", " << track_l1.dxy() << ", " << track_l1.dz() << std::endl;
          std::cerr << "l2 pt, eta, phi, dxy, dz " << track_l2.pt() << ", " << track_l2.eta() << ", " << track_l2.phi()
                    << ", " << track_l2.dxy() << ", " << track_l2.dz()<< std::endl;
        }
        for (const auto& str : { "sv_chi2", "sv_ndof", "sv_prob", "fitted_mass", "fitted_massErr",
                                 "fitted_pt", "vtx_x", "vtx_y", "vtx_z", "vtx_ex", "vtx_ey", "vtx_ez" }) {
          lepton_pair.addUserFloat(str, -1.);
        }
      }
      // cut on the SV info
      if(post_vtx_selection_ && !(*post_vtx_selection_)(lepton_pair) ) continue;
      dilepton_pairs->push_back(lepton_pair);
    }
  }
  evt.put(std::move(dilepton_pairs));
}

using DiMuonBuilder = DiLeptonBuilder<pat::Muon, pat::Muon>;
using DiElectronBuilder = DiLeptonBuilder<pat::Electron, pat::Electron>;
using DiTrackBuilder = DiLeptonBuilder<reco::Track, reco::Track>;
using EleMuBuilder = DiLeptonBuilder<pat::Electron, pat::Muon>;
using MuTrackBuilder = DiLeptonBuilder<pat::Muon, reco::Track>;
using EleTrackBuilder = DiLeptonBuilder<pat::Electron, reco::Track>;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiMuonBuilder);
DEFINE_FWK_MODULE(DiElectronBuilder);
DEFINE_FWK_MODULE(DiTrackBuilder);
DEFINE_FWK_MODULE(EleMuBuilder);
DEFINE_FWK_MODULE(MuTrackBuilder);
DEFINE_FWK_MODULE(EleTrackBuilder);
