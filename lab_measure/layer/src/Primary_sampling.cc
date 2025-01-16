/*
 * Primary_sampling.cc
 *
 *  Created on: Nov 17, 2017
 *      Author: kevin
 */
#include "Primary_sampling.hh"

G4double neutron_energy(){
	G4double x, y;
	while (true){
		x = -log(G4UniformRand());
		y = -log(G4UniformRand());
		if (pow(y - cf_m * (x + 1), 2) <= cf_b * cf_l * x)
			return cf_l * x;
	}
}

G4double gamma_pdf(G4double e){
	G4double y;
	if (e >= 0.085 and e < 0.3)
		y = 38.13 * (e - 0.085) * exp(1.648 * e);
	else if (e >= 0.3 and e < 1.0)
		y = 26.8 * exp(-2.3 * e);
	else if (e >= 1.0 and e <= 8.0)
		y = 8.0 * exp(-1.1 * e);
	else
		G4cerr << "Sampled primary gamma ray energy beyond [0.085, 8] MeV\n";
	return y;
}

G4double gamma_energy(){
	G4double max_pdf = 13.5;
	G4double e;
	while (true){
		// sample an energy
		e = 0.085 +  G4UniformRand() * gamma_range;
		G4double pdf = gamma_pdf(e);
		if (max_pdf * G4UniformRand() <= pdf)
			return e;
	}
}

G4ThreeVector iso_dir(){
	// cosine of polar angle, 0 to 1
	G4double cos_theta = G4UniformRand();
	G4double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

	// azimuthal angle
	G4double phi = CLHEP::twopi * G4UniformRand();

	G4double mu = sin_theta * cos(phi);
	G4double eta = sin_theta * sin(phi);
	G4double xi = cos_theta;

	return G4ThreeVector(mu, eta, xi);



}
