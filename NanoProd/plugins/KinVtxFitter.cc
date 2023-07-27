#include "KinVtxFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"

KinVtxFitter::KinVtxFitter(const std::vector<reco::TransientTrack>& tracks,
                           const std::vector<double>& masses,
                           const std::vector<double>& sigmas):
  n_particles_{masses.size()} {

  KinematicParticleFactoryFromTransientTrack factory;
  std::vector<RefCountedKinematicParticle> particles;

  for(size_t i = 0; i < tracks.size(); ++i) {
    float sigma = sigmas.at(i);
    particles.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_,
        kin_ndof_, sigma
        )
      );
  }

  KinematicParticleVertexFitter kcv_fitter;
  // Move secondary vertex cutoff to the end of the muon system
  kcv_fitter.setTrackerBounds(650., 1000.);
  RefCountedKinematicTree vtx_tree = kcv_fitter.fit(particles);

  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false;
    return;
  }

  vtx_tree->movePointerToTheTop();
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){
    success_ = false;
    return;
  }
  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = vtx_tree->finalStateParticles();
  if(fitted_children_.size() != n_particles_) {
    success_=false;
    return;
  }
  fitted_track_ = fitted_particle_->refittedTransientTrack();
  success_ = true;
}

