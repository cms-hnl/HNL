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

  template<typename Lepton1, typename Lepton2>
  struct DiLeptonBuilderHelper {
    static const Lepton2& getSecondLepton(const edm::Handle<edm::View<Lepton1>>& leptons1,
                                          const edm::Handle<edm::View<Lepton2>>& leptons2,
                                          size_t index, bool use_first)
    {
      if(use_first)
        throw std::runtime_error("use_first = true, but Lepton1 type != Lepton2 type");
      return leptons2->at(index);
    }

    static bool isSameObject(const Lepton1& lepton1, const Lepton2& lepton2) { return false; }
  };

  template<typename Lepton>
  struct DiLeptonBuilderHelper<Lepton, Lepton> {
    static const Lepton& getSecondLepton(const edm::Handle<edm::View<Lepton>>& leptons1,
                                         const edm::Handle<edm::View<Lepton>>& leptons2,
                                         size_t index, bool use_first)
    {
      if(use_first)
        return leptons1->at(index);
      return leptons2->at(index);
    }

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
              consumes<Lepton1Collection>( cfg.getParameter<edm::InputTag>("src")));
    } else {
      src1_ = std::make_unique<edm::EDGetTokenT<Lepton1Collection>>(
              consumes<Lepton1Collection>( cfg.getParameter<edm::InputTag>("src1")));
      src2_ = std::make_unique<edm::EDGetTokenT<Lepton2Collection>>(
              consumes<Lepton2Collection>( cfg.getParameter<edm::InputTag>("src2")));
    }

    if(cfg.exists("srcVeto")) {
      src_veto_ = std::make_unique<edm::EDGetTokenT<edm::View<reco::Candidate>>>(
                  consumes<edm::View<reco::Candidate>>( cfg.getParameter<edm::InputTag>("srcVeto")));
    }

    if(cfg.exists("deltaR_thr")) {
      deltaR_thr_ = cfg.getParameter<double>("deltaR_thr");
    }

    if(cfg.exists("l1l2Interchangeable")) {
      l1l2_interchangeable_ = cfg.getParameter<bool>("l1l2Interchangeable");
      if(l1l2_interchangeable_ && !l1l2HaveSameType)
        throw std::runtime_error("Inconsistent config. l1l2Interchangeable = true, but Lepton1 type != Lepton2 type.");
    }
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
  std::unique_ptr<edm::EDGetTokenT<edm::View<reco::Candidate>>> src_veto_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  bool l1l2_interchangeable_{false};
  double deltaR_thr_{0.3};
};


template<typename Lepton1, typename Lepton2>
void DiLeptonBuilder<Lepton1, Lepton2>::produce(edm::StreamID, edm::Event& evt, edm::EventSetup const& setup) const {

  const TransientTrackBuilder& tt_builder = setup.getData(ttkToken_);

  //input
  edm::Handle<Lepton1Collection> leptons1;
  if(src1_)
    evt.getByToken(*src1_, leptons1);
  edm::Handle<Lepton2Collection> leptons2;
  if(src2_)
    evt.getByToken(*src2_, leptons2);
  edm::Handle<edm::View<reco::Candidate> > veto_leptons;
  if(src_veto_)
    evt.getByToken(*src_veto_, veto_leptons);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

  const size_t n_l1 = leptons1->size();
  const size_t n_l2 = src2_ ? leptons2->size() : n_l1;
  for(size_t l1_idx = 0; l1_idx < n_l1; ++l1_idx) {
    const Lepton1& l1 = leptons1->at(l1_idx);
    if(l1_selection_ && !(*l1_selection_)(l1)) continue;
    if(src_veto_ && dR_match(l1, *veto_leptons, deltaR_thr_)) continue;

    const size_t l2_start = l1l2_interchangeable_ && !src2_ ? l1_idx + 1 : 0;
    for(size_t l2_idx = l2_start; l2_idx < n_l2; ++l2_idx) {
      const Lepton2& l2 = Helper::getSecondLepton(leptons1, leptons2, l2_idx, src2_.get() != nullptr);
      if(Helper::isSameObject(l1, l2)) continue;
      if(l2_selection_ && !(*l2_selection_)(l2)) continue;
      if(src_veto_ && dR_match(l2, *veto_leptons, deltaR_thr_)) continue;

      pat::CompositeCandidate lepton_pair;
      lepton_pair.setP4(getP4(l1) + getP4(l2));
      lepton_pair.setCharge(l1.charge() + l2.charge());
      lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(l1, l2));

      // Put the lepton passing the corresponding selection
      if(!l1l2_interchangeable_ || l1.pt() >= l2.pt()) {
        lepton_pair.addUserInt("l1_idx", l1_idx);
        lepton_pair.addUserInt("l2_idx", l2_idx);
      } else {
        lepton_pair.addUserInt("l1_idx", l2_idx);
        lepton_pair.addUserInt("l2_idx", l1_idx);
      }
      // // Adding user cands would be helpful, but not possible for tracks
      // lepton_pair.addUserCand("l1", l1_ptr );
      // lepton_pair.addUserCand("l2", l2_ptr );

      // before making the SV, cut on the info we have
      if(pre_vtx_selection_ && !(*pre_vtx_selection_)(lepton_pair) ) continue;

      reco::TransientTrack tt_l1(tt_builder.build(getTrack(l1)));
      reco::TransientTrack tt_l2(tt_builder.build(getTrack(l2)));

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
        // std::cout << e.what() << std::endl;
        // std::cout << "l1 pt, eta, phi, dxy, dz " << l1.pt() << ", " << l1.eta() << ", " << l1.phi() << ", " << getTrack(l1).dxy() << ", " << getTrack(l1).dz() << std::endl;
        // std::cout << "l2 pt, eta, phi, dxy, dz " << l2.pt() << ", " << l2.eta() << ", " << l2.phi() << ", " << getTrack(l2).dxy() << ", " << getTrack(l2).dz()<< std::endl;
        for (const auto& str : {"sv_chi2", "sv_ndof", "sv_prob", "fitted_mass", "fitted_massErr", "fitted_pt", "vtx_x", "vtx_y", "vtx_z", "vtx_ex", "vtx_ey", "vtx_ez"}) {
          lepton_pair.addUserFloat(str, -1.);
        }
      }

      // cut on the SV info
      if(post_vtx_selection_ && !(*post_vtx_selection_)(lepton_pair) ) continue;
      ret_value->push_back(lepton_pair);
    }
  }

  evt.put(std::move(ret_value));
}

using DiMuonBuilder = DiLeptonBuilder<pat::Muon, pat::Muon>;
using DiElectronBuilder = DiLeptonBuilder<pat::Electron, pat::Electron>;
using DiTrackBuilder = DiLeptonBuilder<reco::Track, reco::Track>;
using MuEleBuilder = DiLeptonBuilder<pat::Muon, pat::Electron>;
using MuTrackBuilder = DiLeptonBuilder<pat::Muon, reco::Track>;
using EleTrackBuilder = DiLeptonBuilder<pat::Electron, reco::Track>;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiMuonBuilder);
DEFINE_FWK_MODULE(DiElectronBuilder);
DEFINE_FWK_MODULE(DiTrackBuilder);
DEFINE_FWK_MODULE(MuEleBuilder);
DEFINE_FWK_MODULE(MuTrackBuilder);
DEFINE_FWK_MODULE(EleTrackBuilder);
