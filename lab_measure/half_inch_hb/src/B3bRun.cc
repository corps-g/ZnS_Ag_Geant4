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
// $Id: B3bRun.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B3bRun.cc
/// \brief Implementation of the B3bRun class

#include "B3bRun.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3bRun::B3bRun()
 : G4Run(), 
   totalEnergy(0.),
   totalEnergy2(0.),
   protonEnergy(0.0),
   protonEnergy2(0.0)
{
	G4SDManager* SDM = G4SDManager::GetSDMpointer();
	totalEnergyID = SDM -> GetCollectionID("ZnSAg/edepTotal");
	protonEnergyID = SDM -> GetCollectionID("ZnSAg/edepProton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3bRun::~B3bRun()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3bRun::RecordEvent(const G4Event* event)
{
  
  //Hits collections
  //  
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(!HCE) return;
   
  G4THitsMap<G4double>* evtMapProton =
    static_cast<G4THitsMap<G4double>*>(HCE->GetHC(protonEnergyID));
  //G4cout << "Event " << evtNb << ", size of cryst event map = " << 
	//        evtMap->entries() << G4endl;
               
  std::map<G4int,G4double*>::iterator itr;
  for (itr = evtMapProton->GetMap()->begin();
		  itr != evtMapProton->GetMap()->end(); itr++) {
	  G4double tmp = *(itr -> second);
	  protonEnergy += tmp;
	  protonEnergy2 += tmp * tmp;
  }
  G4THitsMap<G4double>* evtMapTotal =
		  static_cast<G4THitsMap<G4double>*>(HCE->GetHC(totalEnergyID));
  for (itr = evtMapTotal -> GetMap() -> begin();
		  itr != evtMapTotal -> GetMap()->end(); itr++) {
	  G4double tmp = *(itr -> second);
	  totalEnergy += tmp;
	  totalEnergy2 += tmp * tmp;
  }  
  G4Run::RecordEvent(event);      
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3bRun::Merge(const G4Run* aRun)
{
  const B3bRun* localRun = static_cast<const B3bRun*>(aRun);
  totalEnergy += localRun -> totalEnergy;
  totalEnergy2 += localRun -> totalEnergy2;
  protonEnergy += localRun -> protonEnergy;
  protonEnergy2 += localRun -> protonEnergy2;
  G4Run::Merge(aRun); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
