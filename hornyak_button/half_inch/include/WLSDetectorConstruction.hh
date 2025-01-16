#ifndef WLSDetectorConstruction_h
#define WLSDetectorConstruction_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4RotationMatrix.hh"
#include "G4GenericMessenger.hh"

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

#include <vector>
typedef std::vector<G4double> vec_double;
typedef std::vector<vec_double> vec2_double;
typedef std::vector<vec2_double> vec3_double;

typedef std::vector<G4int> vec_int;
typedef std::vector<vec_int> vec2_int;

class WLSDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    WLSDetectorConstruction(G4double, G4double);
    virtual ~WLSDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();

    virtual void ConstructSDandField();

    G4Material* FindMaterial(G4String);

    // sample the center of a grain
    vec_double sample(G4int x, G4int y, G4int z, G4int numx, G4int numy,
    		G4double dx, G4double dy, G4double startx, G4double starty,
			G4double startz);

    // compute the distance^2 of two grains
    G4double distance2(vec_double, vec_double);

    // compare the sampled center with the neighbor grains
    G4bool compare(vec_double, vec2_double);

    // find the neighbor grain(s) need to check overlap
    vec2_double find_neighbors(G4int x, G4int y, G4int z,
    		vec_double sc, vec3_double currentz, vec3_double lastz,
			G4double startx, G4double starty, G4double startz);

    G4double compute_startx(G4int, G4double);
    G4double compute_starty(G4int, G4double);
    G4double compute_startz(G4int);

  private:

    G4VPhysicalVolume* fPhysiWorld;
    G4VPhysicalVolume* fPhysiPMMA1;
    G4VPhysicalVolume* fPhysiPMMA2;
    G4VPhysicalVolume* scinPV;
    G4VPhysicalVolume* pmtPV;

    G4LogicalVolume* aZSALV;
    G4LogicalVolume* aPMMALV1;
    G4LogicalVolume* aPMMALV2;
    G4LogicalVolume* scinLV;
    G4LogicalVolume* pmtLV;

    WLSMaterials* fMaterials;
    G4double lenx, leny, lenz, diameter, rZSA, dZSA, wr, rhoZSA, rhoPMMA, inch;

    G4double absLen, sigmaalpha;
    G4double dZSA2;
    vec2_int zlayer;  // grain distribution info

    G4int nx, nx1, ny, ny1, nz;
    G4double dx0, dx1, dy0, dy1, dz;

    G4double edgez, edgey, edgex;
    G4double startx0, starty0, startz0;
};

#endif
