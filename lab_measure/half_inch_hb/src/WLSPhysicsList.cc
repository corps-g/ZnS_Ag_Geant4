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
// $Id: WLSPhysicsList.cc 78066 2013-12-03 11:08:36Z gcosmo $
//
/// \file optical/wls/src/WLSPhysicsList.cc
/// \brief Implementation of the WLSPhysicsList class
//
//
#include "WLSPhysicsList.hh"
#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//#include "G4PhysListFactory.hh"
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT_HP.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"


#include "G4ProcessTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

#include "G4RadioactiveDecayPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4OpticalPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhysicsList::WLSPhysicsList(G4String physName) : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();
//    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = NULL;
    if (physName == "QGSP_BERT_HP") {
       phys = new QGSP_BERT_HP;
    } else if (physName == "FTFP_BERT_HP"){
       phys = new FTFP_BERT_HP;
    }
//    if (factory.IsReferencePhysList(physName)) {
//       phys = factory.GetReferencePhysList(physName);
//       if(!phys)G4Exception("WLSPhysicsList::WLSPhysicsList","InvalidSetup",
//                            FatalException,"PhysicsList does not exist");
//    }
    for (G4int i = 0; ; ++i) {
       G4VPhysicsConstructor* elem =
                  const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
       if (elem == NULL) break;
       G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
       RegisterPhysics(elem);
    }
	G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
//	opticalPhysics->SetScintillationByParticleType(true);

  RegisterPhysics( opticalPhysics );

    fAbsorptionOn = true;
    
    //This looks complex, but it is not:
    //Get from base-class the pointer of the phsyicsVector
    //to be used. Remember: G4VModularPhysicsList is now a split class.
    //Why G4VModularPhysicsList::RegisterPhysics method is not used instead?
    //If possible we can remove this...
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhysicsList::~WLSPhysicsList()
{}

void WLSPhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  //
  SetCutsWithDefault();

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


