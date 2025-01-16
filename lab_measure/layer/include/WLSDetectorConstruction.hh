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
// $Id: WLSDetectorConstruction.hh 77487 2013-11-25 10:15:04Z gcosmo $
//
/// \file optical/wls/include/WLSDetectorConstruction.hh
/// \brief Definition of the WLSDetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSDetectorConstruction_h
#define WLSDetectorConstruction_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4RotationMatrix.hh"

class G4Box;
class G4Tubs;
class G4EllipticalTube;

class G4LogicalVolume;
class G4VPhysicalVolume;

class WLSMaterials;
class G4Material;

class WLSDetectorMessenger;

class WLSPhotonDetSD;

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"

class WLSDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    WLSDetectorConstruction(G4double, G4double, G4int);
    virtual ~WLSDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    virtual void ConstructSDandField();

    G4Material* FindMaterial(G4String);

  private:
    void DefineCommands();

    G4VPhysicalVolume* fPhysiWorld;
    G4VPhysicalVolume* fPhysiPMMA;
    G4VPhysicalVolume* fPhysizsa;

    G4LogicalVolume* aZSALV;
    G4LogicalVolume* aPMMALV;

    WLSMaterials* fMaterials;
    // lenx, leny: lengths of total button
    // lenPMMA, lenZSA: length of one PMMA and ZSA layer
    G4double lenx, leny, lenPMMA, lenZSA;
    G4int nUnit;
};

#endif
