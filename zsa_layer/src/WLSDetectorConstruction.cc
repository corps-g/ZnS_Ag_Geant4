//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: WLSDetectorConstruction.cc 84718 2014-10-20 07:40:45Z gcosmo $
//
/// \file optical/wls/src/WLSDetectorConstruction.cc
/// \brief Implementation of the WLSDetectorConstruction class
//
//
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSDetectorConstruction::WLSDetectorConstruction(G4double len1,
		G4double len2, G4int n)
: fPhysiWorld(NULL)
{
	lenx = 0.1255 * 2.0 * cm;
	leny = 0.4445 * 2.0 * cm;
	lenPMMA = len1;
	lenZSA = len2;
	nUnit = n;
	G4cout << "lenpmma = " << lenPMMA / cm << " cm, lenZSA = " << lenZSA / um
			<< " um, number of units = " << nUnit << G4endl;
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

G4VPhysicalVolume* WLSDetectorConstruction::ConstructDetector()
{
	G4double lenUnit = lenPMMA + lenZSA;
	G4double lenz = nUnit * (lenPMMA + lenZSA);
	//--------------------------------------------------
	// World
	//--------------------------------------------------
	G4VSolid* solidWorld = new G4Box("World", lenx, leny, lenz+.1*cm);

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
	// PMMA layer
	//--------------------------------------------------
	G4VSolid* solidPMMA = new G4Box("solidPMMA", .5 * lenx, .5 * leny,
			                        .5 * lenPMMA);
	aPMMALV = new G4LogicalVolume(solidPMMA, FindMaterial("G4_PLEXIGLASS"),
			                      "pmmaLV");

	//--------------------------------------------------
	// ZSA layer
	//--------------------------------------------------
	G4VSolid* solidzsa = new G4Box("zsaSolid", 0.5 * lenx, 0.5 *leny,
			                        0.5 * lenZSA);
	aZSALV = new G4LogicalVolume(solidzsa, FindMaterial("ZnSAg"),
			                      "ZnS(AG)LV");

	//--------------------------------------------------
	// Place layers
	//--------------------------------------------------
	for(G4int i=0; i<nUnit; i++){
		G4double zPMMA = .5 * lenPMMA + i * lenUnit;
		G4double zZSA = zPMMA + .5 * lenUnit;
		new G4PVPlacement(0, G4ThreeVector(0., 0., zPMMA),
				          aPMMALV,
						  "pmmaPV" + std::to_string(i),
						  fLogicWorld,
						  false,
						  0,
						  false);
		new G4PVPlacement(0, G4ThreeVector(0., 0., zZSA),
				          aZSALV,
						  "zsaPV" + std::to_string(i),
						  fLogicWorld,
						  false,
						  0,
						  false);
	}

	//--------------------------------------------------
	//Optical surface property
	//--------------------------------------------------
	G4double opE[] = {1. * eV, 2. * eV, 3. * eV, 4. * eV};
	G4double reflectivity1[] = {1., 1., 1., 1.};
	G4int n = 4;

	//Optical skin surface surrounding ZnS(Ag) and PMMA layers
	G4OpticalSurface* opSurface = new G4OpticalSurface("polish",
			                         unified,
									 polished,
									 dielectric_dielectric);
	G4OpticalSurface* groundsurface = new G4OpticalSurface("ground",
				                         unified,
										 ground,
										 dielectric_dielectric);

	G4MaterialPropertiesTable* mptRel1 = new G4MaterialPropertiesTable();
	mptRel1 -> AddProperty("REFLECTIVITY", opE, reflectivity1, n);
	opSurface ->SetMaterialPropertiesTable(mptRel1);
	groundsurface ->SetMaterialPropertiesTable(mptRel1);

	new G4LogicalSkinSurface("ZnS(Ag)SkinSur", aZSALV, groundsurface);
	new G4LogicalSkinSurface("PMMASkinSurface", aPMMALV, opSurface);
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

	// Deposited energy by alpha
	G4VPrimitiveScorer* edepAlpha = new G4PSEnergyDeposit("edepAlpha");
	G4SDParticleFilter* alphaFilter = new	G4SDParticleFilter("alphaFilter");
	alphaFilter ->	add("alpha");
	edepAlpha -> SetFilter(alphaFilter);
	fZnSAgE -> RegisterPrimitive(edepAlpha);

	// E by generic ion
	G4VPrimitiveScorer* edepIon = new G4PSEnergyDeposit("edepIon");
	G4SDParticleFilter* ionFilter = new G4SDParticleFilter("ionFilter");
	ionFilter -> add("GenericIon");
	edepIon -> SetFilter(ionFilter);
	fZnSAgE -> RegisterPrimitive(edepIon);

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
	G4Material* material = G4Material::GetMaterial(name,true);
	return material;
}
