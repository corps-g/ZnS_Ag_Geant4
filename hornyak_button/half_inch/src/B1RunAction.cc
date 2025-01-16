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
// $Id: B1RunAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParameterManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <math.h>
#include "B3bRun.hh"
#include "Analysis.hh"
#include <string>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction()
{ 
	// Analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	std::vector<G4double> edges = {0., 0.9};
	for (int i = 0; i < 400; i++)
		edges.push_back((i+1) * 5. + .1);
	analysisManager->CreateH1("0","Total number of detected op",
			edges);
	analysisManager->SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction(){
	delete G4AnalysisManager::Instance();
}

G4Run* B1RunAction::GenerateRun()
{ return new B3bRun; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
	// inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);

}

vecG4double B1RunAction::cal_ave(G4double var, G4double var2, G4int n){
	vecG4double ans(2);
	G4double ave = var / n, ave2 = var2 / n;
	ans[0] = ave;
	G4double std;
	std = sqrt((ave2 - ave * ave) / (n - 1));
	ans[1] = std / ave;
	return ans;
}

void B1RunAction::EndOfRunAction(const G4Run* run)
{
	// Analysis manager
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	if (IsMaster()) {
		G4cout
		<< G4endl
		<< "--------------------End of Global Run-----------------------\n";
		G4int nofEvents = run->GetNumberOfEvent();
		if (nofEvents == 0) return;
		G4cout << "Mean number of optical photon = " <<
				analysisManager -> GetH1(0) -> mean() << G4endl;

		G4cout << "rms = " <<
				analysisManager -> GetH1(0) -> rms() << G4endl;

		const B3bRun* b3Run = static_cast<const B3bRun*>(run);
		G4double totalEnergy = b3Run -> GetTotalDepoE();
		G4double totalEnergy2 = b3Run -> GetTotalDepoE2();
		G4double protonEnergy = b3Run -> GetProtonDepoE();
		G4double protonEnergy2 = b3Run -> GetProtonDepoE2();
		vecG4double total, proton;
		total = cal_ave(totalEnergy, totalEnergy2, nofEvents);
		proton = cal_ave(protonEnergy, protonEnergy2, nofEvents);

		G4cout
		<< "Ave. total deposited energy in ZnS(Ag), relative error = "
		<< total[0] << " MeV,  " << total[1] << "\n"
		<< "Ave. proton deposited energy in ZnS(Ag), relative error = "
		<< proton[0] << " MeV, " << proton[1] << G4endl;
	}
	else {
		G4cout
		<< G4endl
		<< "--------------------End of Local Run------------------------\n";
		G4cout << "Mean number of op = "
					<< analysisManager -> GetH1(0) -> mean() << G4endl;
		G4cout << "rms = " << analysisManager -> GetH1(0) -> rms() << G4endl;
	}
	analysisManager -> Write();
	analysisManager -> CloseFile();
}
