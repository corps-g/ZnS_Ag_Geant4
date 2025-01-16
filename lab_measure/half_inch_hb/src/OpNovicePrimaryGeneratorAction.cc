#include "OpNovicePrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

OpNovicePrimaryGeneratorAction::OpNovicePrimaryGeneratorAction()
:G4VUserPrimaryGeneratorAction(),
 gun(0),
 neutron(0),
 gamma(0)
{
	gun = new G4ParticleGun(1);
	auto ptable = G4ParticleTable::GetParticleTable();
	neutron = ptable -> FindParticle("neutron");
	gamma = ptable -> FindParticle("gamma");

	// point souce at fixed position
	// half-inch Hornyak button sample
	G4double lenz = 0.5 * 2.54 * cm;
	// source 2-cm away from Hornyak button surface, origin of the coordinate is
	// at the center of the Hornyak button
	G4double srcZ = -0.5 * lenz - 2.0 * cm;
	gun -> SetParticlePosition(G4ThreeVector(0.0, 0.0, srcZ));
}


OpNovicePrimaryGeneratorAction::~OpNovicePrimaryGeneratorAction()
{
	if (gun)
		delete gun;
}


void OpNovicePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4int ng;
	// sample number of gamma rays in 2 or 3 for 2.13 gamma rays per neutron
	if (G4UniformRand() < 0.13)
		ng = 3;
	else
		ng = 2;

	// one neutron per event
	gun -> SetParticleDefinition(neutron);
	gun -> SetParticleEnergy(neutron_energy());
	gun -> SetParticleMomentumDirection(iso_dir());
	gun -> GeneratePrimaryVertex(anEvent);

	// 2 or 3 gamma rays per event
	for (G4int i = 0; i < ng; i++){
		gun -> SetParticleDefinition(gamma);
		gun -> SetParticleEnergy(gamma_energy());
		gun -> SetParticleMomentumDirection(iso_dir());
		gun -> GeneratePrimaryVertex(anEvent);
	}
}
