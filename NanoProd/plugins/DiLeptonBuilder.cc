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


template<typename Lepton>
class DiLeptonBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<Lepton> LeptonCollection;

  explicit DiLeptonBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<LeptonCollection>( cfg.getParameter<edm::InputTag>("src") )} {
       produces<pat::CompositeCandidateCollection>();
    }

  ~DiLeptonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<Lepton> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton> l2_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<LeptonCollection> src_;
};

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

template<typename Lepton>
void DiLeptonBuilder<Lepton>::produce(edm::StreamID, edm::Event& evt, edm::EventSetup const& setup) const {
  
  edm::ESHandle<TransientTrackBuilder> tt_builder ;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", tt_builder);


  //input
  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(src_, leptons);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  
  for(size_t l1_idx = 0; l1_idx < leptons->size(); ++l1_idx) {
    edm::Ptr<Lepton> l1_ptr(leptons, l1_idx);
    if(!l1_selection_(*l1_ptr)) continue; 
    
    for(size_t l2_idx = l1_idx + 1; l2_idx < leptons->size(); ++l2_idx) {
      edm::Ptr<Lepton> l2_ptr(leptons, l2_idx);
      if(!l2_selection_(*l2_ptr)) continue;

      pat::CompositeCandidate lepton_pair;
      lepton_pair.setP4(getP4(*l1_ptr) + getP4(*l2_ptr));
      lepton_pair.setCharge(l1_ptr->charge() + l2_ptr->charge());
      lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));
      // Put the lepton passing the corresponding selection
      lepton_pair.addUserInt("l1_idx", l1_idx );
      lepton_pair.addUserInt("l2_idx", l2_idx );
      // // Adding user cands would be helpful, but not possible for tracks
      // lepton_pair.addUserCand("l1", l1_ptr );
      // lepton_pair.addUserCand("l2", l2_ptr );
      if( !pre_vtx_selection_(lepton_pair) ) continue; // before making the SV, cut on the info we have


      reco::TransientTrack tt_l1(tt_builder->build(getTrack(*l1_ptr)));
      reco::TransientTrack tt_l2(tt_builder->build(getTrack(*l2_ptr)));

      KinVtxFitter fitter(
        {tt_l1, tt_l2},
        {LEP_MASS, LEP_MASS},
        {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
        );

      lepton_pair.setVertex(
		     reco::Candidate::Point(
					    fitter.fitted_vtx().x(),
					    fitter.fitted_vtx().y(),
					    fitter.fitted_vtx().z()
					    )
		     );

      lepton_pair.addUserFloat("sv_chi2", fitter.chi2());
      lepton_pair.addUserFloat("sv_ndof", fitter.dof()); // float??
      lepton_pair.addUserFloat("sv_prob", fitter.prob());
      lepton_pair.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1.);
      lepton_pair.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1.);
      lepton_pair.addUserFloat("vtx_x", lepton_pair.vx());
      lepton_pair.addUserFloat("vtx_y", lepton_pair.vy());
      lepton_pair.addUserFloat("vtx_z", lepton_pair.vz());
      lepton_pair.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      lepton_pair.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      lepton_pair.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));
     
 // if needed, add here more stuff

      // cut on the SV info
      if( !post_vtx_selection_(lepton_pair) ) continue;
      ret_value->push_back(lepton_pair);
    }
  }
  
  evt.put(std::move(ret_value));
}

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
typedef DiLeptonBuilder<pat::Muon> DiMuonBuilder;
typedef DiLeptonBuilder<pat::Electron> DiElectronBuilder;
typedef DiLeptonBuilder<reco::Track> DiTrackBuilder;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiMuonBuilder);
DEFINE_FWK_MODULE(DiElectronBuilder);
DEFINE_FWK_MODULE(DiTrackBuilder);


