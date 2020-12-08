// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "PhysicsTools/PatUtils/interface/MiniIsolation.h"


class TrackIsoTableProducer : public edm::global::EDProducer<> {
public:
  explicit TrackIsoTableProducer(const edm::ParameterSet& iConfig)
      : name_(iConfig.getParameter<std::string>("name")),
        src_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dsa"))),
        srcPF_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pf"))) {
    produces<nanoaod::FlatTable>();
  }

  ~TrackIsoTableProducer() override{};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("dsa", edm::InputTag("selectedDSAMuons"))->setComment("input dsa muon collection");
    desc.add<edm::InputTag>("pf", edm::InputTag("packedPFCandidates"))
        ->setComment("input pf candidates");
    desc.add<std::string>("name", "DSAMuons")->setComment("name of the dsa muon nanoaod::FlatTable we are extending with isolation");
    descriptions.add("dsaTable", desc);
  }

private:
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  std::string name_;
  edm::EDGetTokenT<std::vector<reco::Track>> src_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> srcPF_;
};


namespace {
   
   pat::PFIsolation getMiniPFIsolation(const pat::PackedCandidateCollection& pfcands,
                                  float eta,
                                  float phi,
                                  float drcut,
                                  float ptthresh=0.5,
                                  float deadcone_ch=0.0001,
                                  float deadcone_pu=0.01,
                                  float deadcone_ph=0.01,
                                  float deadcone_nh=0.01,
                                  float dZ_cut=0.) {
     float chiso = 0, nhiso = 0, phiso = 0, puiso = 0;
     for (auto const &pc : pfcands) {
       float dr2 = deltaR2(eta, phi, pc.eta(), pc.phi());
       if (dr2 > drcut * drcut)
         continue;
       int id = pc.pdgId();
       float pt = pc.pt();
       if (std::abs(id) == 211) {
         bool fromPV = (pc.fromPV() > 1 || fabs(pc.dz()) < dZ_cut);
         if (fromPV && dr2 > deadcone_ch * deadcone_ch) {
           // if charged hadron and from primary vertex, add to charged hadron isolation
           chiso += pt;
         } else if (!fromPV && pt > ptthresh && dr2 > deadcone_pu * deadcone_pu) {
           // if charged hadron and NOT from primary vertex, add to pileup isolation
           puiso += pt;
         }
       }
       // if neutral hadron, add to neutral hadron isolation
       if (std::abs(id) == 130 && pt > ptthresh && dr2 > deadcone_nh * deadcone_nh)
         nhiso += pt;
       // if photon, add to photon isolation
       if (std::abs(id) == 22 && pt > ptthresh && dr2 > deadcone_ph * deadcone_ph)
         phiso += pt;
     }
 
     return pat::PFIsolation(chiso, nhiso, phiso, puiso);
   }
}

// ------------ method called to produce the data  ------------
void TrackIsoTableProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<std::vector<reco::Track>> muons;
  iEvent.getByToken(src_, muons);
  edm::Handle<std::vector<pat::PackedCandidate>> pfs;
  iEvent.getByToken(srcPF_, pfs);

  int precision = 10;

  unsigned int ncand = muons->size();

  std::vector<float> ch_iso(ncand, 0.), nh_iso(ncand, 0.), ph_iso(ncand, 0.), pu_iso(ncand, 0.);
  
  for (unsigned int i = 0; i < ncand; ++i) {
    const reco::Track& mu = (*muons)[i];

    pat::PFIsolation iso = getMiniPFIsolation(*pfs, mu.eta(), mu.phi(), 0.4);

    ch_iso[i] = iso.chargedHadronIso();
    nh_iso[i] = iso.neutralHadronIso();
    ph_iso[i] = iso.photonIso();
    pu_iso[i] = iso.puChargedHadronIso();
  }

  auto tab = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);

  tab->addColumn<float>(
      "chIso",
      ch_iso,
      "charged isolation pT sum",
      nanoaod::FlatTable::FloatColumn,
      precision);
  tab->addColumn<float>(
      "nhIso",
      nh_iso,
      "neutral isolation pT sum",
      nanoaod::FlatTable::FloatColumn,
      precision);
  tab->addColumn<float>(
      "phIso",
      ph_iso,
      "photon isolation pT sum",
      nanoaod::FlatTable::FloatColumn,
      precision);
  tab->addColumn<float>(
      "puIso",
      pu_iso,
      "pileup isolation pT sum",
      nanoaod::FlatTable::FloatColumn,
      precision);  
  iEvent.put(std::move(tab));
}


#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(TrackIsoTableProducer);
