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


WLSDetectorConstruction::WLSDetectorConstruction(G4double lz, G4double awr)
: fPhysiWorld(NULL){
	lenx = 0.1255 * 2.0 * cm;
	leny = 0.4445 * 2.0 * cm;
	lenz = lz;
	rZSA = 20.0 * um;
	dZSA = 2.0 * rZSA;
	dZSA2 = dZSA * dZSA;
	wr = awr;
	rhoZSA = 4.09;
	rhoPMMA = 1.19;

	// edge of the scintillation volume
	edgez = -0.5 * lenz;
	edgey = -0.5 * leny;
	edgex = -0.5 * lenx;

	// start corordinates of the ZnS(Ag) grain
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
	fMaterials = WLSMaterials::GetInstance();
	return ConstructDetector();
}

vec_double WLSDetectorConstruction::sample(G4int x, G4int y, G4int z,
		G4int numx, G4int numy,
		G4double dx, G4double dy,
		G4double startx,
		G4double starty,
		G4double startz){
	G4double cx, cy, cz;

	if (x == 0 or x == numx - 1){
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
	G4double dist2;
	dist2 = pow((g1[0] - g2[0]), 2) + pow((g1[1] - g2[1]), 2) +
			pow((g1[2] - g2[2]), 2);
	return dist2;
}

G4bool WLSDetectorConstruction::compare(vec_double sc, vec2_double neighbors){
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
	// x, y, z: the sampled grain cell indexes
	// sc: the sampled center
	// currentz: placed grain centers in current z layer
	// lastz: placed grain centers in z-1 layer
	// start*: start corordinate of this cell in * axis

	//  centers of the neighbor grains need to be check
	vec2_double neighbors;

	// the current z layer

	if (x != 0){
		// if x is not the first cell, the left grain may needs to be check
		if (sc[0] - startx < dZSA)
			neighbors.push_back(currentz[y][x - 1]);
	}

	if (y != 0){
		// if y is not the bottom layer, the grains in y-1 layer may need to be
		// check
		if (sc[1] - starty < dZSA){
			// cell width in x axis of y-1 layer
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
			// left and right edges of the sampled center
			G4double leftx = corx - dZSA, rightx = corx + dZSA;
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

	// the last z layer
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

			// the sample center edges along y axis in z-1 layer
			G4int upj = G4int((cory + dZSA - edgey) / dy);
			G4int downj = G4int((cory - dZSA - edgey) / dy);

			// global boundary check
			if (upj >= numy)
				upj = numy - 1;
			if (downj < 0)
				downj = 0;

			for(G4int j = downj; j <= upj; j++){
				// determine dx in this y layer
				if (zlayer[z - 1][j]){
					// nx1
					dx = dx1;
					numx = nx1;
				}
				else{
					dx = dx0;
					numx = nx;
				}

				// determine grains in x axis in this y layer
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
	G4cout << "lenz = " << lenz / cm << " cm\n";
	G4cout << "zsa wr = " << wr << "\n";

	G4VSolid* solidWorld = new G4Box("World", lenx, leny, lenz);

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
	// pmma
	//--------------------------------------------------
	G4VSolid* solidPMMA = new G4Box("PMMASolid", 0.5 * lenx, 0.5 *leny, 0.5 * lenz);
	aPMMALV = new G4LogicalVolume(solidPMMA, FindMaterial("G4_PLEXIGLASS"), "pmmaLV");
	fPhysiPMMA = new G4PVPlacement(0, G4ThreeVector(),
			aPMMALV, "pmmaPV", fLogicWorld, false, 0, true);

	// solid and logical volume of a ZnS(Ag) grain
	G4VSolid* solidZSA = new G4Orb("solidZSA", rZSA);
	aZSALV = new G4LogicalVolume(solidZSA, FindMaterial("ZnSAg"), "ZnS(AG)LV");

	//--------------------------------------------------
	// Place ZnS(Ag) grains
	//--------------------------------------------------

	// Calculate distribution data of ZnS(Ag) grains
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
	G4double nZSA = volZSA / (4. / 3. * 3.141592653589793 * rZSA * rZSA * rZSA);

	 /* The nZSA grains are distributed uniformly in the scintillation volume.
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
		The choose of length axis (z) is because it is the longest dimension,
		which has the largest number of layers.
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
								aPMMALV,
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
					// place a dummy grain to indicate it is a void cell
					currentz[y][x] = dummy;
				}
			} // end of x loop
		} // end of y loop
		lastz = currentz;
	} // end of z loop

	G4cout << "\nFinish ZnS(Ag) grain placement" << G4endl;

	// check the real weight ratio
	G4double relative_error = (nZSA - cnt) / nZSA;
	G4cout << "\nRelative error of weight ratio = " << relative_error << G4endl;
	G4cout << "The sampled weight ratio = " << cnt / nZSA * wr << G4endl;

	//--------------------------------------------------
	//Optical surface property
	//--------------------------------------------------
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

	// ground surface
	G4OpticalSurface* groundsurface = new G4OpticalSurface("genericOpS",
			unified,
			ground,
			dielectric_dielectric);
	groundsurface -> SetMaterialPropertiesTable(mptRel1);

	new G4LogicalSkinSurface("ground_zsa_surface", aZSALV, groundsurface);
	new G4LogicalSkinSurface("polish_pmma_surface", aPMMALV, opSurface);

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

	// E by electron
	G4VPrimitiveScorer* edepElec = new G4PSEnergyDeposit("edepElec");
	G4SDParticleFilter* elecFilter = new G4SDParticleFilter("elecFilter");
	elecFilter -> add("e-");
	edepElec -> SetFilter(elecFilter);
	fZnSAgE -> RegisterPrimitive(edepElec);

	// Total deposited energy
	G4VPrimitiveScorer* edepTotal = new G4PSEnergyDeposit("edepTotal");
	fZnSAgE -> RegisterPrimitive(edepTotal);
	SetSensitiveDetector("ZnS(AG)LV", fZnSAgE);
}

G4Material* WLSDetectorConstruction::FindMaterial(G4String name) {
	G4Material* material = G4Material::GetMaterial(name, true);
	return material;
}

