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
// $Id: WLSSteppingAction.cc 75292 2013-10-30 09:25:15Z gcosmo $
//
/// \file optical/wls/src/WLSSteppingAction.cc
/// \brief Implementation of the WLSSteppingAction class
//
//
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"

#include "WLSSteppingAction.hh"
#include "WLSDetectorConstruction.hh"

#include "G4ParticleTypes.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"

#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>
#include "B1EventAction.hh"
#include "G4SolidStore.hh"
#include "G4Box.hh"
#include "G4VSolid.hh"

WLSSteppingAction::WLSSteppingAction(B1EventAction* eventAction, G4double lz)
:fEventAction(eventAction){
	halfz = lz * .5;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSSteppingAction::~WLSSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
	G4Track* theTrack = theStep->GetTrack();
	// Need to know if this is an optical photon
	if(theTrack->GetDefinition()
			!= G4OpticalPhoton::OpticalPhotonDefinition()) return;

	G4StepPoint* thePrePoint  = theStep->GetPreStepPoint();
	G4StepPoint* thePostPoint = theStep->GetPostStepPoint();

	G4VPhysicalVolume* thePrePV  = thePrePoint->GetPhysicalVolume();
	G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

	G4String thePrePVname  = " ";
	G4String thePostPVname = " ";

	if (thePostPV) {
		thePrePVname  = thePrePV->GetName();
		thePostPVname = thePostPV->GetName();
	}

	// Retrieve the status of the photon
	G4OpBoundaryProcessStatus theStatus = Undefined;
	G4ProcessManager* OpManager =
			G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

	if (OpManager) {
		G4int MAXofPostStepLoops =
				OpManager->GetPostStepProcessVector()->entries();
		G4ProcessVector* fPostStepDoItVector =
				OpManager->GetPostStepProcessVector(typeDoIt);

		for ( G4int i=0; i<MAXofPostStepLoops; i++) {
			G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
			fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
			if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
		}
	}

	switch (theStatus) {
	case FresnelRefraction:
		if (thePostPVname == "World") {
			G4double locz = thePostPoint -> GetPosition().z();
			if (locz > - halfz){
				// Not emitted from input surface
				fEventAction -> increase_total();}
			theTrack->SetTrackStatus(fStopAndKill);
			return;
		}
		break;
	default: break;
	}
}
