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
// $Id: B1ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B1ActionInitialization.cc
/// \brief Implementation of the B1ActionInitialization class

#include "B1ActionInitialization.hh"
#include "B1RunAction.hh"
#include "B1EventAction.hh"
#include "OpNovicePrimaryGeneratorAction.hh"
#include "WLSSteppingAction.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ActionInitialization::B1ActionInitialization(G4double len1, G4double len2,
		G4int n)
 : G4VUserActionInitialization()
{
	lenPMMA = len1;
	lenZSA = len2;
	nUnit = n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ActionInitialization::~B1ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ActionInitialization::BuildForMaster() const
{
  B1RunAction* runAction = new B1RunAction(nUnit);
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ActionInitialization::Build() const
{
  SetUserAction(new OpNovicePrimaryGeneratorAction());

  B1RunAction* runAction = new B1RunAction(nUnit);
  SetUserAction(runAction);
  
  B1EventAction* eventAction = new B1EventAction(runAction, nUnit);
  SetUserAction(eventAction);
  
  SetUserAction(new WLSSteppingAction(eventAction, lenPMMA, lenZSA, nUnit));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
