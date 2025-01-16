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
// $Id: WLSSteppingAction.hh 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/include/WLSSteppingAction.hh
/// \brief Definition of the WLSSteppingAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSSteppingAction_h
#define WLSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "B1EventAction.hh"


class G4Track;
class G4StepPoint;

class G4OpBoundaryProcess;

class WLSSteppingAction : public G4UserSteppingAction
{
  public:

    WLSSteppingAction(B1EventAction* eventAction, G4double lz);
    virtual ~WLSSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
 
  private:

    // number of photons that reach the end
    B1EventAction*  fEventAction;
    G4OpBoundaryProcess* fOpProcess;
    G4double halfz;
};

#endif
