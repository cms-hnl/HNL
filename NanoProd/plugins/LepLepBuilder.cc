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
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "KinVtxFitter.h"

constexpr float LEP_SIGMA = 0.0000001;
constexpr float LEP_MASS = 0.1056584;


template<typename Lepton1, typename Lepton2>
class DiLeptonBuilder : public edm::global::EDProducer<> {

public:
  typedef std::vector<Lepton1> Lepton1Collection;
  typedef std::vector<Lepton2> Lepton2Collection;

  explicit DiLeptonBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src1_{consumes<Lepton1Collection>( cfg.getParameter<edm::InputTag>("src1") )},
    src2_{consumes<Lepton2Collection>( cfg.getParameter<edm::InputTag>("src2") )} {
       produces<pat::CompositeCandidateCollection>();
    }

  ~DiLeptonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<Lepton1> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton2> l2_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<Lepton1Collection> src1_;
  const edm::EDGetTokenT<Lepton2Collection> src2_;
};

namespace {
    template<typename Lepton>
    const reco::Track& getTrack(const Lepton& lep) {
    return *lep.bestTrack();
    }

    template<>
    const reco::Track& getTrack<reco::Track>(const reco::Track& lep) {
    return lep;
    }

    template<typename Lepton>
    const math::XYZTLorentzVector getP4(const Lepton& lep) {
    return lep.p4();
    }

    template<>
    const math::XYZTLorentzVector getP4<reco::Track>(const reco::Track& lep) {
    float px = lep.px(), py = lep.py(), pz = lep.pz();
    return math::XYZTLorentzVector{px, py, pz, sqrt(px*px + py*py + pz*pz)};
    }
}
template<typename Lepton1, typename Lepton2>
void DiLeptonBuilder<Lepton1, Lepton2>::produce(edm::StreamID, edm::Event& evt, edm::EventSetup const& setup) const {
  
  edm::ESHandle<TransientTrackBuilder> tt_builder ;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", tt_builder);

  //input
  edm::Handle<Lepton1Collection> leptons1;
  evt.getByToken(src1_, leptons1);
  edm::Handle<Lepton2Collection> leptons2;
  evt.getByToken(src2_, leptons2);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

  for(size_t l1_idx = 0; l1_idx < leptons1->size(); ++l1_idx) {
    const Lepton1& l1 = leptons1->at(l1_idx);
    if(!l1_selection_(l1)) continue; 
    
    for(size_t l2_idx = 0; l2_idx < leptons2->size(); ++l2_idx) {
      const Lepton2& l2 = leptons2->at(l2_idx);
      if(!l2_selection_(l2)) continue;

      pat::CompositeCandidate lepton_pair;
      lepton_pair.setP4(getP4(l1) + getP4(l2));
      lepton_pair.setCharge(l1.charge() + l2.charge());
      lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(l1, l2));
      // Put the lepton passing the corresponding selection
      lepton_pair.addUserInt("l1_idx", l1_idx );
      lepton_pair.addUserInt("l2_idx", l2_idx );
      // // Adding user cands would be helpful, but not possible for tracks
      // lepton_pair.addUserCand("l1", l1_ptr );
      // lepton_pair.addUserCand("l2", l2_ptr );
      if( !pre_vtx_selection_(lepton_pair) ) continue; // before making the SV, cut on the info we have

      reco::TransientTrack tt_l1(tt_builder->build(getTrack(l1)));
      reco::TransientTrack tt_l2(tt_builder->build(getTrack(l2)));

      try {
        KinVtxFitter fitter(
          {tt_l1, tt_l2},
          {LEP_MASS, LEP_MASS},
          {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
          );
        if (fitter.success()) {
          lepton_pair.setVertex(
            reco::Candidate::Point(
                  fitter.fitted_vtx().x(),
                  fitter.fitted_vtx().y(),
                  fitter.fitted_vtx().z()
                  )
            );
        } else {
          lepton_pair.setVertex(reco::Candidate::Point());
        }

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
        std::cout << e.what() << std::endl;
        std::cout << "l1 pt, eta, phi, dxy, dz " << l1.pt() << ", " << l1.eta() << ", " << l1.phi() << ", " << getTrack(l1).dxy() << ", " << getTrack(l1).dz() << std::endl;
        std::cout << "l2 pt, eta, phi, dxy, dz " << l2.pt() << ", " << l2.eta() << ", " << l2.phi() << ", " << getTrack(l2).dxy() << ", " << getTrack(l2).dz()<< std::endl;
        for (const auto& str : {"sv_chi2", "sv_ndof", "sv_prob", "fitted_mass", "fitted_massErr", "fitted_pt", "vtx_x", "vtx_y", "vtx_z", "vtx_ex", "vtx_ey", "vtx_ez"}) {
          lepton_pair.addUserFloat(str, -1.);
        } 
      } 

      // cut on the SV info
      if( !post_vtx_selection_(lepton_pair) ) continue;
      ret_value->push_back(lepton_pair);
    }
  }
  
  evt.put(std::move(ret_value));
}

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
typedef DiLeptonBuilder<pat::Muon, pat::Electron> MuEleBuilder;
typedef DiLeptonBuilder<pat::Muon, reco::Track> MuTrackBuilder;
typedef DiLeptonBuilder<pat::Electron, reco::Track> EleTrackBuilder;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuEleBuilder);
DEFINE_FWK_MODULE(MuTrackBuilder);
DEFINE_FWK_MODULE(EleTrackBuilder);


