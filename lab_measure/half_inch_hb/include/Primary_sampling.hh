/*
 * Primary_sampling.hh
 *
 *  Created on: Nov 17, 2017
 *      Author: kevin
 */

#ifndef HALF_INCH_HB_INCLUDE_PRIMARY_SAMPLING_HH_
#define HALF_INCH_HB_INCLUDE_PRIMARY_SAMPLING_HH_
#include "Randomize.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <cmath>
#include "G4SystemOfUnits.hh"
/*
 * Sample the energy and direction of the primary particles from Cf252
 * spontaneous fission, i.e., fission neutron and fission gamma ray.
 */

// constants for cf252 watt fission spectrum, MCNP manual
static constexpr G4double cf_a = 1.025;
static constexpr G4double cf_b = 2.926;

static constexpr G4double cf_k = 1.0 + (cf_b / 8.0 / cf_a);
static constexpr G4double cf_l = (cf_k + sqrt(cf_k * cf_k - 1)) / cf_a;
static constexpr G4double cf_m = cf_a * cf_l - 1.0;


// Sample neutron energy from Cf252 Watt fission spectrum
// https://wci.llnl.gov/codes/tart/media/pdf/UCRL-TR-203351.pdf
G4double neutron_energy();


// Sample the fission gamma ray energy
// literature: Simulation of Neutron and Gamma Ray Emission from Fission
//             and Photofission

// Compute the gamma energy pdf, i.e., a piecewise function based on sampled
// gamma energy
G4double gamma_pdf(G4double);

static G4double gamma_range = 8.0 - 0.085;
// Based on the fission gamma PDF, sample the gamma energy by rejection method
G4double gamma_energy();


// Sample half isotropic (+z) direction vector
G4ThreeVector iso_dir();



#endif /* HALF_INCH_HB_INCLUDE_PRIMARY_SAMPLING_HH_ */
