#include "G4ios.hh"
#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"

#include "WLSMaterials.hh"
#include "WLSDetectorConstruction.hh"

#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"
#include <math.h>
#include <assert.h>
#include <string>
#include "G4Orb.hh"

#include "Randomize.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"


WLSDetectorConstruction::WLSDetectorConstruction(G4double len, G4double sm)
:fPhysiWorld(NULL)
{
	inch = 2.54 * cm;
	lenx = 5. / 8. * inch;
	leny = 7. / 64. * inch;
	lenz = 1. * inch;

	// global boundary coordinates of the scintillation volume
	edgez = -0.5 * lenz;
	edgey = -0.5 * leny;
	edgex = -0.5 * lenx;

	diameter = .75 / 2. / cos(10. / 180. * pi) * 2.54 * cm * 2.;
	rZSA = 20.0 * um;
	dZSA = 2.0 * rZSA;
	dZSA2 = dZSA * dZSA;
	wr = 0.05;
	rhoZSA = 4.09;
	rhoPMMA = 1.19;
	absLen = len;
	sigmaalpha = sm;

	// start coordinates for grain position sampling, do not cross the global
	// boundary
	startx0 = edgex + rZSA;
	starty0 = edgey + rZSA;
	startz0 = edgez + rZSA;
}

WLSDetectorConstruction::~WLSDetectorConstruction() {
	if (fMaterials) delete fMaterials;
}

G4VPhysicalVolume* WLSDetectorConstruction::Construct()
{
	if (fPhysiWorld) {
		G4GeometryManager::GetInstance()->OpenGeometry();
		G4PhysicalVolumeStore::GetInstance()->Clean();
		G4LogicalVolumeStore::GetInstance()->Clean();
		G4SolidStore::GetInstance()->Clean();
		G4LogicalSkinSurface::CleanSurfaceTable();
		G4LogicalBorderSurface::CleanSurfaceTable();
	}
	fMaterials = WLSMaterials::GetInstance(absLen);
	return ConstructDetector();
}

vec_double WLSDetectorConstruction::sample(G4int x, G4int y, G4int z,
		G4int numx, G4int numy,
		G4double dx, G4double dy,
		G4double startx,
		G4double starty,
		G4double startz){
	G4double cx, cy, cz;

	/* Sample the center position of a ZnS(Ag) grain.
	 * Parameters:
	 * 	x, y, z: the cell indexes in which the grain sit
	 * 	numx, numy: number of y cells in this z layer and number of x cells in
	 * 	            the y&z layer (number of z divisions is constant, nz)
	 * 	dx, dy: lengths of this cell in x and y axis
	 * 	startx, starty, startz: start coordinates of this cell
	 */
	if (x == 0 or x == numx - 1){
		// if the cell is the first or the last one of this x layer, reduced the
		// sampling distance to avoid crossing the global boundary
		// if it is the first cell, the start coordinate has already moved
		// into the scintillation volume from the global boundary by rZSA
		cx = startx + G4UniformRand() * (dx - rZSA);
	}
	else{
		cx = startx + G4UniformRand() * dx;
	}

	if (y == 0 or y == numy - 1){
		cy = starty + G4UniformRand() * (dy - rZSA);
	}
	else{
		cy = starty + G4UniformRand() * dy;
	}

	if (z == 0 or z == nz - 1){
		cz = startz + G4UniformRand() * (dz - rZSA);
	}
	else{
		cz = startz + G4UniformRand() * dz;
	}
	vec_double sc = {cx, cy, cz};
	return sc;
}

G4double WLSDetectorConstruction::distance2(vec_double g1, vec_double g2){
	/*
	 * compute the distance^2 between two grains
	 */
	G4double dist2;
	dist2 = pow((g1[0] - g2[0]), 2) + pow((g1[1] - g2[1]), 2) +
			pow((g1[2] - g2[2]), 2);
	return dist2;
}

G4bool WLSDetectorConstruction::compare(vec_double sc, vec2_double neighbors){
	/*
	 * Compare the distance between the sampled grain with the neighbor grains
	 * that need to check.
	 * If no conflict, return true, otherwise, return false.
	 */
	G4double dis2;
	for(G4int n = 0; n < neighbors.size(); n++){
		if (neighbors[n][0] < edgex)
			// a void cell
			continue;
		dis2 = distance2(sc, neighbors[n]);
		if (dis2 < dZSA2)
			return false;
	}
	return true;
}

vec2_double WLSDetectorConstruction::find_neighbors(G4int x, G4int y, G4int z,
		vec_double sc, vec3_double currentz, vec3_double lastz,
		G4double startx, G4double starty, G4double startz){
	/*
	 * For a sampled grain, find the neighbor placed grains that need to check.
     * Parameters:
		   x, y, z: the sampled grain cell indexes
		   sc: the sampled grain center
		   currentz: placed grain centers in current z layer
		   lastz: placed grain centers in the z-1 layer
		   start*: start coordinate of this cell in * axis
	 * Return: a vector containing the centers of the placed grains that need to
	 *         be checked.
	 */

	//  centers of the neighbor grains need to be check
	vec2_double neighbors;

	// placed grains in the current z layer
	if (x != 0){
		// if x is not the first cell, the left grain may needs to be check
		// The [x-1] grain needs to be checked if the distance between the x
		// coordinate of the sampled grain and the cell x edge < dZSA.
		// Otherwise, the two grains can not conflict.
		if (sc[0] - startx < dZSA)
			neighbors.push_back(currentz[y][x - 1]);
	}

	// the placed grains in the y-1 layer of current z layer
	if (y != 0){
		// if y is not the bottom layer, the grains in y-1 layer may need to be
		// check
		if (sc[1] - starty < dZSA){
			// cell width in x axis of y-1 layer, which has two candidate values
			G4double dx;
			G4int numx;
			if (zlayer[z][y-1]){
				dx = dx1;
				numx = nx1;
			}
			else{
				dx = dx0;
				numx = nx;
			}

			G4double corx = sc[0];
			// left and right edges of the sampled center, the sensitive region
			G4double leftx = corx - dZSA, rightx = corx + dZSA;
			// locate the placed grains in the y-1 layer need to be check
			G4int li = G4int((leftx - edgex) / dx);
			G4int ri = G4int((rightx - edgex) / dx);

			// global boundary check
			if (li < 0)
				li = 0;
			if (ri >= numx)
				ri = numx - 1;

			for (G4int i = li; i <= ri; i++){
				neighbors.push_back(currentz[y - 1][i]);
			}
		}
	}

	// the grains in the z-1 layer
	if (z != 0){
		// z!=0 layers, need to consider z-1 layer
		if (sc[2] - startz < dZSA){
			G4double corx = sc[0], cory = sc[1];

			// xy layer information of z-1 layer
			G4double dx, dy;
			G4int numx, numy;
			numy = zlayer[z-1].size();
			if (numy == ny)
				dy = dy0;
			else
				dy = dy1;

			// the sample center edges along y axis in z-1 layer, the sensitive
			// region
			G4int upj = G4int((cory + dZSA - edgey) / dy);
			G4int downj = G4int((cory - dZSA - edgey) / dy);

			// global boundary check
			if (upj >= numy)
				upj = numy - 1;
			if (downj < 0)
				downj = 0;

			for(G4int j = downj; j <= upj; j++){
				// determine dx in this z, y layer
				if (zlayer[z - 1][j]){
					// nx1
					dx = dx1;
					numx = nx1;
				}
				else{
					dx = dx0;
					numx = nx;
				}

				// determine grains in x axis in this z, y layer
				G4int li = G4int((corx - dZSA - edgex) / dx);
				G4int ri = G4int((corx + dZSA - edgex) / dx);

				// global boundary check
				if (li < 0)
					li = 0;
				if (ri >= numx)
					ri = numx - 1;

				for(G4int i = li; i <= ri; i++)
					neighbors.push_back(lastz[j][i]);
			}
		}
	}
	return neighbors;
}

G4double WLSDetectorConstruction::compute_startx(G4int x, G4double dx){
	// Compute the start x coordinate of a cell for sampling.
	// For the first cell, move the startx inside the sctintillation volume by
	// rZSA, i.e., startx0.
	if (x == 0)
		return startx0;
	else
		return edgex + x * dx;
}

G4double WLSDetectorConstruction::compute_starty(G4int y, G4double dy){
	if (y == 0)
		return starty0;
	else
		return edgey + y * dy;
}

G4double WLSDetectorConstruction::compute_startz(G4int z){
	if (z == 0)
		return startz0;
	else
		return edgez + z * dz;
}

G4VPhysicalVolume* WLSDetectorConstruction::ConstructDetector()
{
	//--------------------------------------------------
	// World
	//--------------------------------------------------
	G4VSolid* solidWorld = new G4Box("World", diameter, diameter, lenz);

	G4LogicalVolume* fLogicWorld = new G4LogicalVolume(solidWorld,
			FindMaterial("G4_AIR"),
			"World");

	fPhysiWorld = new G4PVPlacement(0,
			G4ThreeVector(),
			fLogicWorld,
			"World",
			0,
			false,
			0,
			true);

	//--------------------------------------------------
	// Button
	//--------------------------------------------------

	// light guide 1, part of a half cylinder
	G4VSolid* tmpsolidPMMA1 = new G4Tubs("tmplight_guide1", 0.0, 0.5 * diameter,
			                        0.5 * lenz, 0., pi);
	G4double subh = tan(10. / 180. * pi) * .75 * .5 * 2.54 * cm;
	G4VSolid* subbox = new G4Box("subbox", diameter, subh, .6 * lenz);
	G4SubtractionSolid* solidPMMA1 = new G4SubtractionSolid("light_guide1",
			tmpsolidPMMA1, subbox);
	aPMMALV1 = new G4LogicalVolume(solidPMMA1, FindMaterial("G4_PLEXIGLASS"),
			                      "pmmaLV1");
	fPhysiPMMA1 = new G4PVPlacement(0, G4ThreeVector(0., .5 * leny - subh, 0.),
									   aPMMALV1, "pmmaPV1", fLogicWorld,
									   false, 0, true);

	// light guide 2, part of a half cylinder
	G4VSolid* tmpsolidPMMA2 = new G4Tubs("tmplight_guide2", 0.0, 0.5 * diameter,
				                        0.5 * lenz, pi, pi);
	G4SubtractionSolid* solidPMMA2 = new G4SubtractionSolid("light_guide1",
				tmpsolidPMMA2, subbox);
	aPMMALV2 = new G4LogicalVolume(solidPMMA2, FindMaterial("G4_PLEXIGLASS"),
					                      "pmmaLV2");
	fPhysiPMMA2 = new G4PVPlacement(0, G4ThreeVector(0., -(0.5 * leny - subh), 0.),
				                       aPMMALV2, "pmmaPV2", fLogicWorld,
									   false, 0, true);

	// scintillation volume
	G4VSolid* scin = new G4Box("scin", .5 * lenx, .5 * leny, .5 * lenz);
	scinLV = new G4LogicalVolume(scin, FindMaterial("G4_PLEXIGLASS"), "scinLV");
	scinPV = new G4PVPlacement(0, G4ThreeVector(),
			                   scinLV, "scinPV", fLogicWorld, false, 0, true);

	// dummy PMT, for tally purpose
	G4VSolid* pmt = new G4Tubs("pmt", 0., 1.5 * diameter / 2.,
			                    0.5 * 0.1 * lenz, 0., 2. * pi);
	pmtLV = new G4LogicalVolume(pmt, FindMaterial("G4_AIR"), "pmtLV");
	pmtPV = new G4PVPlacement(0,
			                  G4ThreeVector(0., 0.,
			                		  0.5 * lenz + 0.1 * 0.5 * lenz),
							  pmtLV, "pmtPV", fLogicWorld, false, 0, true);

	// source position sampling
	// union of the scintillation volume and two light guides
	// the source particles are born in a big plane, z = srcz, and the particles
	// are confine to this volume
	G4VSolid* srcScin = new G4Box("srcScin", .5*lenx, .5*leny, .5 * .1 * lenz);

	G4VSolid* tmpsrcPMMA1 = new G4Tubs("tmpsrcPMMA1", 0., .5 * diameter,
			                        .5 * .1 * lenz, 0, pi);
	G4SubtractionSolid* srcPMMA1 = new G4SubtractionSolid("srcPMMA1",
			tmpsrcPMMA1, subbox);

	G4VSolid* tmpsrcPMMA2 = new G4Tubs("tmpsrcPMMA2", 0., .5 * diameter,
			                         .5 * .1 * lenz, pi, pi);
	G4SubtractionSolid* srcPMMA2 = new G4SubtractionSolid("srcPMMA2",
				tmpsrcPMMA2, subbox);
	// displace a little bit to avoid sharing surface
	G4UnionSolid* union1 = new G4UnionSolid("union1", srcScin, srcPMMA1, 0,
			                                G4ThreeVector(0, .5*leny-1.e-6-subh, 0));
	G4UnionSolid* union2 = new G4UnionSolid("union2", union1, srcPMMA2, 0,
			                                G4ThreeVector(0, -(.5*leny-subh) + 1.e-6, 0));

	G4LogicalVolume* srcLV = new G4LogicalVolume(union2,
				                                 FindMaterial("G4_AIR"),
				                                 "srcLV");
	new G4PVPlacement(0, G4ThreeVector(0., 0., -.5*lenz - .05 * lenz),
			          srcLV, "srcPV", fLogicWorld, false, 0, true);


	// solid and logical volume of a ZnS(Ag) grain
	G4VSolid* solidZSA = new G4Orb("solidZSA", rZSA);
	aZSALV = new G4LogicalVolume(solidZSA, FindMaterial("ZnSAg"), "ZnS(AG)LV");

	//--------------------------------------------------
	// Place ZnS(Ag) grains
	//--------------------------------------------------
	// volume of the scintillation volume
	G4double vol = lenx * leny * lenz;

	/* The volume of the ZnS(Ag) grains can be computed based on the 5% weight
	 * ratio, densities of ZnS(Ag) and PMMA, and the modeling of ZnS(Ag) grain
	 * as sphere with radius of 20 um.
	 * wr = (rho_zsa * V_zsa) / (rho_zsa * V_zsa + rho_pmma * V_pmma)   (1)
	 * V_zsa + V_pmma = V_total = lenx * leny * lenz                    (2)
	 * The two equations lead to
	 * V_zsa = wr * rho_pmma * V_total / [(1 - wr) * rho_zsa + wr * rho_pmma]
	*/
	G4double volZSA = wr * rhoPMMA * vol / ((1.0 - wr) * rhoZSA + wr * rhoPMMA);

	// total number of ZnS(Ag) grains in the scintillation volume
	// = V_zsa / V_sphere
	G4double nZSA = volZSA / (4. / 3. * pi * rZSA * rZSA * rZSA);

	/*
	 * The nZSA grains are distributed uniformly in the scintillation volume.
	 * Assume
	 * nz / lenz = nx / lenx = ny / leny = k,                            (3)
	 * where n* is the number of ZnS(Ag) grains along * axis.
	 * nx * ny * nz = nZSA                                               (4)
	 * lenx * leny * lenz = V_total                                      (5)
	 * Eqs. (3), (4), and (5) lead to
	 * k = (nZSA / V_total) ^ (1/3)
	 */

	G4double k = pow(nZSA / vol, 1. / 3.);
	G4double fnz = k * lenz;

	/*
	 * The number of layers in the z axis is an int, nz.
	In each z layer, the number of y layers is sampled to fny.
	In each z, y layer, the number of x layers is sampled to fnx.
	The choose of length axis (z) is because it is the longest dimension, which has
	the largest number of layers.
	Two advantages:
		1) allows enough layer number to sample fny and fnx
		2) effect of the round of (fnz -> nz) is smallest.
   	   	   loss_n = (fnz - nz) * fnx * fny
   	   	   The lost grains are added into xy payer anyway.
	 */
	// constant number of z layers
	nz = G4int(fnz);
	dz = lenz / nz;

	// ensure the layer is enough for a grain
	assert (dz > dZSA);

	// number of grains in a z layer
	G4double fxy = nZSA / nz;

	/*
	 * Distribute the nxy into x and y proportional to the lenx, leny
		fnx / lenx = fny / leny = kxy
		kxy^2 = fnx * fny / (lenx * leny) = fxy / (lenx * leny)
		fnx = kxy * lenx
		fny = kxy * leny
	 */
	G4double kxy = pow(fxy / (lenx * leny), 0.5);
	G4double fny = kxy * leny;
	G4double fnx = kxy * lenx;

	G4cout << G4endl << "fnx = " << fnx << "; fny = " << fny << G4endl;

	// Precompute the parameters in x axis
	nx = G4int(fnx);
	nx1 = nx + 1;
	// probability of number = nx1 grains in this x layer
	G4double diffx = fnx - nx;
	// cell length, evenly divided into nx1 cells, leave a void cell at head and
	// tail repeatedly if there are nx grains
	dx0 = lenx / nx;
	dx1 = lenx / nx1;
	assert (dx1 > dZSA);

	// precompute the parameters in y axis
	ny = G4int(fny);
	ny1 = ny + 1;
	G4double diffy = fny - ny;
	dy0 = leny / ny;
	dy1 = leny / ny1;
	assert (dy1 > dZSA);

	// precompute the grain distribution
	G4int numy, numx;
	for(G4int z = 0; z < nz; z++){
		// sample number of y layers
		if (G4UniformRand() < diffy)
			numy = ny1;
		else
			numy = ny;
		vec_int xylayer(numy, 0);
		for (G4int y = 0; y < numy; y++){
			if (G4UniformRand() < diffx)
				// nx1 grains in this z, y layer
				xylayer[y] = 1;
			// default is 0, i.e., nx
		}
		zlayer.push_back(xylayer);
	}

	// dummy grain center to indicate it is a void
	vec_double dummy = {-lenx, 0.0, 0.0};

	// track the grains in the current z layer and the last z layer for overlap
	// check
	vec3_double currentz, lastz;

	// variables used in the placement
	G4String name;
	G4int cnt = 0, voidcnt = 0;
	G4int fail_time, max_time = 100;
	G4double dx, dy;
	// center of the sampled grain
	vec_double sc;

	G4double startx, starty, startz;
	// place the grains
	for(G4int z = 0; z < nz; z++){
		// compute startz
		startz = compute_startz(z);

		// y layer information
		vec_int xylayer = zlayer[z];
		numy = xylayer.size();
		// y layer width
		if (numy == ny)
			dy = dy0;
		else
			dy = dy1;

		// set the number of xy grains in this z layer
		currentz.resize(numy);
		for (G4int i = 0; i < numy; i++){
			if (xylayer[i])
				currentz[i].resize(nx1);
			else
				currentz[i].resize(nx);
		}

		for (G4int y = 0; y < numy; y++){
			// compute starty
			starty = compute_starty(y, dy);

			// number of grains in this z, y layer
			if (xylayer[y]){
				numx = nx1;
				dx = dx1;
			}
			else{
				numx = nx;
				dx = dx0;
			}

			for (G4int x = 0; x < numx; x++){
				// compute startx
				startx = compute_startx(x, dx);

				fail_time = 0;
				while(fail_time < max_time){
					// sample a center
					sc = sample(x, y, z, numx, numy, dx, dy,
							startx, starty, startz);

					// find neighbors
					vec2_double neighbors = find_neighbors(x, y, z, sc, currentz,
							lastz, startx, starty, startz);

					if (compare(sc, neighbors)){
						name = std::to_string(z) +
								std::to_string(y) +
								std::to_string(x);
						new G4PVPlacement(0,
								G4ThreeVector(sc[0], sc[1], sc[2]),
								aZSALV,
								name,
								scinLV,
								false,
								0);
						currentz[y][x] = sc;
						cnt ++;
						break;
					}
					fail_time ++;
				}
				if (fail_time == max_time){
					voidcnt ++;
					currentz[y][x] = dummy;
				}
			} // end of x loop
		} // end of y loop
		lastz = currentz;
	} // end of z loop

	G4cout << "\nFinish ZnS(Ag) grain placement" << G4endl;

	// check the real weight ratio
	G4double relative_error = (nZSA - cnt) / nZSA;
	G4cout << "Relative error of weight ratio = " << relative_error << G4endl;
	G4cout << "The sampled weight ratio = " << cnt / nZSA * wr << G4endl;

	/*--------------------------------------------------
	  Optical surface property. Information can be found at P206 at Geant4
	  manual (version 10.3) and the literature:
	    Identifying key surface parameters for optical photon transport in
        GEANT4/GATE simulations
	  --------------------------------------------------*/
	G4double opE[] = {1. * eV, 2. * eV, 3. * eV, 4. * eV};
	G4double reflectivity1[] = {1., 1., 1., 1.};
	G4int n = 4;

	// polish surface
	G4OpticalSurface* opSurface = new G4OpticalSurface("genericOpS",
			                         unified,
									 polished,
									 dielectric_dielectric);

	G4MaterialPropertiesTable* mptRel1 = new G4MaterialPropertiesTable();
	mptRel1 -> AddProperty("REFLECTIVITY", opE, reflectivity1, n);
	opSurface ->SetMaterialPropertiesTable(mptRel1);

	// white paint surface
	G4OpticalSurface* paintsurface = new G4OpticalSurface("genericOpS",
                                                          unified,
														  polishedfrontpainted,
			                                              dielectric_dielectric);
	paintsurface -> SetMaterialPropertiesTable(mptRel1);

	// ground surface
	G4OpticalSurface* groundsurface = new G4OpticalSurface("genericOpS",
				                         unified,
										 ground,
										 dielectric_dielectric, sigmaalpha);
	groundsurface -> SetMaterialPropertiesTable(mptRel1);

	// ground skin surface of the ZnS(Ag) grains
	new G4LogicalSkinSurface("ZnS(Ag)SkinSur", aZSALV, groundsurface);

	// Scintillation volume and light guides have white paint surfaces (default).
	new G4LogicalSkinSurface("pmmaSkinSur", aPMMALV1, paintsurface);
	new G4LogicalSkinSurface("pmmaSkinSur", aPMMALV2, paintsurface);
	new G4LogicalSkinSurface("scin_outer_surface", scinLV, paintsurface);

	// border surface, which takes precedence over the skin surface (Geant4
	// manual). Overwrite the white paint skin surfaces at the couple surfaces.

	// Coupled surfaces between light guides and scintillation volume are
	// defined as polished ones.
	new G4LogicalBorderSurface("scin_to_light_guide1", scinPV, fPhysiPMMA1,
			                   opSurface);
	new G4LogicalBorderSurface("light_guide1_to_scin", fPhysiPMMA1, scinPV,
			opSurface);
	new G4LogicalBorderSurface("scin_to_light_guide2", scinPV, fPhysiPMMA2,
			opSurface);
	new G4LogicalBorderSurface("light_guide2_to_scin", fPhysiPMMA2, scinPV,
			opSurface);

	// Coupled surfaces between (scintillation volume + light guides) to PMT
	// are polished surfaces.
	new G4LogicalBorderSurface("light_guide1_to_pmt", fPhysiPMMA1, pmtPV,
			                   opSurface);
	new G4LogicalBorderSurface("light_guide2_to_pmt", fPhysiPMMA2, pmtPV,
				               opSurface);
	new G4LogicalBorderSurface("scin_to_pmt", scinPV, pmtPV,
					           opSurface);
	return fPhysiWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::ConstructSDandField() {
	G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
	G4MultiFunctionalDetector* fZnSAgE = new G4MultiFunctionalDetector("ZnSAg");
	// Deposited energy by proton
	G4VPrimitiveScorer* edepProton = new G4PSEnergyDeposit("edepProton");
	G4SDParticleFilter* protonFilter = new	G4SDParticleFilter("protonFilter");
	protonFilter ->	add("proton");
	edepProton -> SetFilter(protonFilter);
	fZnSAgE -> RegisterPrimitive(edepProton);
	// Total deposited energy
	G4VPrimitiveScorer* edepTotal = new G4PSEnergyDeposit("edepTotal");
	fZnSAgE -> RegisterPrimitive(edepTotal);
	SetSensitiveDetector("ZnS(AG)LV", fZnSAgE);
}

G4Material* WLSDetectorConstruction::FindMaterial(G4String name) {
	G4Material* material = G4Material::GetMaterial(name,true);
	return material;
}

