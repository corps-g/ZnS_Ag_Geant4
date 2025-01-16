#ifndef OpNovicePrimaryGeneratorAction_h
#define OpNovicePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "Primary_sampling.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

/*
 * Model a isotropic point Cf252 neutron and gamma source. There are one neutron
 *  and 2.13 gamma rays with specific fission spectra. The number of gamma rays
 *  were sampled to be 2 or 3 to preserve 2.13.
 *     Position: 2-cm from Hornyak button surface (z), x, y at the center.
 * Liter: http://www.sciencedirect.com/science/article/pii/S0168900215003654
 */

class OpNovicePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    OpNovicePrimaryGeneratorAction();
    virtual ~OpNovicePrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* gun;
    G4ParticleDefinition* neutron;
    G4ParticleDefinition* gamma;
};
#endif
