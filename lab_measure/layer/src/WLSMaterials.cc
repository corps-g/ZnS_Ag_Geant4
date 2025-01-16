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
// $Id: WLSMaterials.cc 82854 2014-07-14 09:08:25Z gcosmo $
//
/// \file optical/wls/src/WLSMaterials.cc
/// \brief Implementation of the WLSMaterials class
//
//
#include "WLSMaterials.hh"

#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include <vector>
#include "G4Element.hh"

WLSMaterials* WLSMaterials::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSMaterials::WLSMaterials()
{
  fNistMan = G4NistManager::Instance();

  fNistMan->SetVerbose(2);

  CreateMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSMaterials::~WLSMaterials()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSMaterials* WLSMaterials::GetInstance()
{
  if (fInstance == 0)
    {
      fInstance = new WLSMaterials();
    }
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* WLSMaterials::GetMaterial(const G4String material)
{
  G4Material* mat =  fNistMan->FindOrBuildMaterial(material);

  if (!mat) mat = G4Material::GetMaterial(material);
  if (!mat) {
     std::ostringstream o;
     o << "Material " << material << " not found!";
     G4Exception("WLSMaterials::GetMaterial","",
                 FatalException,o.str().c_str());
  }

  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSMaterials::CreateMaterials()
{
  // Materials Definitions
  // =====================

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  fAir = fNistMan->FindOrBuildMaterial("G4_AIR");

  //--------------------------------------------------
  // PMMA
  //--------------------------------------------------
  fPMMA = fNistMan -> FindOrBuildMaterial("G4_PLEXIGLASS");

  //--------------------------------------------------
  // Polystyrene
  //--------------------------------------------------
  fPolystyrene = fNistMan -> FindOrBuildMaterial("G4_POLYSTYRENE");

  //--------------------------------------------------
  // ZnS(Ag), 2 wt% Ag
  //--------------------------------------------------

  ZnS = new G4Material("ZnS", 4.09 * g/cm3, 2);
  zinc = fNistMan -> FindOrBuildElement(30);
  sulfide = fNistMan -> FindOrBuildElement(16);
  ZnS -> AddElement(zinc, 1);
  ZnS -> AddElement(sulfide, 1);
  silver = fNistMan -> FindOrBuildElement(47);
  fZnSAg = new G4Material("ZnSAg", 4.09 * g/cm3, 2);
  fZnSAg -> AddMaterial(ZnS, 98.0 * perCent);
  fZnSAg -> AddElement(silver, 2.0 * perCent);
 
  //--------------------------------------------------
  // Aluminium
  //--------------------------------------------------

  //fNistMan->FindOrBuildMaterial("G4_Al");

  //--------------------------------------------------------
  // Generate & Add Material Properties Table
  //--------------------------------------------------------

  G4int n;

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  G4double photonEnergy[] =
    {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
     2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
     2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
     2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
     2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
     2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
     2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
     3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
     3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
     3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

  G4double refractiveIndex[] =
  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};

  assert(sizeof(refractiveIndex) == sizeof(photonEnergy));
  n = sizeof(refractiveIndex) / sizeof(G4double);

  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, n);

  fAir->SetMaterialPropertiesTable(mpt);

  //--------------------------------------------------
  //  PMMA
  //--------------------------------------------------

  G4double ope[] = {1.14482167528 * eV, 1.15226939994 * eV, 1.15981466261 * eV,
		  1.16636112355 * eV, 1.17409268403 * eV, 1.18192743025 * eV,
		  1.18986744178 * eV, 1.19675856596 * eV, 1.20489978069 * eV,
		  1.21315251891 * eV, 1.22151908801 * eV, 1.23000185946 * eV,
		  1.23736714005 * eV, 1.24632275264 * eV, 1.25490068252 * eV,
		  1.26346873976 * eV, 1.27228514554 * eV, 1.28122545658 * eV,
		  1.29029230339 * eV, 1.2994883915 * eV, 1.3088165041 * eV,
		  1.31813935183 * eV, 1.32773813914 * eV, 1.33747775009 * eV,
		  1.3473613066 * eV, 1.35739202357 * eV, 1.36757321237 * eV,
		  1.3777551665 * eV, 1.38824529653 * eV, 1.39889639437 * eV,
		  1.40971219367 * eV, 1.42069654444 * eV, 1.43185341764 * eV,
		  1.44301894126 * eV, 1.45453058931 * eV, 1.46622738213 * eV,
		  1.47811382252 * eV, 1.49019456049 * eV, 1.50247439933 * eV,
		  1.51477321238 * eV, 1.52746319371 * eV, 1.54036759142 * eV,
		  1.55349188614 * eV, 1.56684174691 * eV, 1.5804230393 * eV,
		  1.59403686594 * eV, 1.60809581625 * eV, 1.6224049651 * eV,
		  1.6369710514 * eV, 1.65180105826 * eV, 1.66667814805 * eV,
		  1.6820538249 * eV, 1.69771583504 * eV, 1.71367225201 * eV,
		  1.72993145574 * eV, 1.74650214725 * eV, 1.76314259717 * eV,
		  1.78035880863 * eV, 1.79791455094 * eV, 1.81581996826 * eV,
		  1.83408561292 * eV, 1.85272246613 * eV, 1.87145943295 * eV,
		  1.89086758324 * eV, 1.91068250012 * eV, 1.93091710689 * eV,
		  1.95158488011 * eV, 1.9726998796 * eV, 1.99395605393 * eV,
		  2.01600304769 * eV, 2.03854303573 * eV, 2.06159274082 * eV,
		  2.08516965074 * eV, 2.10929206249 * eV, 2.13361189869 * eV,
		  2.15887493354 * eV, 2.18474339089 * eV, 2.2112392979 * eV,
		  2.23838576337 * eV, 2.26620704502 * eV, 2.29430398655 * eV,
		  2.3235417435 * eV, 2.35353430966 * eV, 2.38431129679 * eV,
		  2.41590388607 * eV, 2.44834493351 * eV, 2.48117245213 * eV,
		  2.51540246365 * eV, 2.55059015497 * eV, 2.58677628694 * eV,
		  2.62400396684 * eV, 2.66231881969 * eV, 2.7011805541 * eV,
		  2.74179981055 * eV, 2.78365934964 * eV, 2.82681685894 * eV,
		  2.87133365987 * eV, 2.91727499843 * eV, 2.96400161208 * eV,
		  3.01298146861 * eV, 3.06360730005 * eV};
  G4double ri[] = {1.48146079143, 1.4815180213, 1.48157597525,
		  1.48162624845, 1.481685625, 1.48174581086,
		  1.48180683729, 1.48185983901, 1.48192251235,
		  1.48198612032, 1.48205069763, 1.48211627998,
		  1.48217332118, 1.48224281448, 1.48230952864,
		  1.4823763258, 1.48244523733, 1.4825153134,
		  1.48258659433, 1.48265912175, 1.48273293865,
		  1.48280697444, 1.48288348451, 1.48296142129,
		  1.48304083396, 1.48312177344, 1.48320429248,
		  1.4832871961, 1.48337301502, 1.4834605828,
		  1.48354996044, 1.48364121127, 1.48373440108,
		  1.48382818339, 1.48392542786, 1.48402482355,
		  1.48412644734, 1.48423037926, 1.48433670264,
		  1.48444388591, 1.48455521796, 1.48466921182,
		  1.48478596589, 1.48490558292, 1.48502817019,
		  1.48515196882, 1.48528078996, 1.48541293126,
		  1.4855485208, 1.48568769266, 1.4858284582,
		  1.48597516465, 1.48612589287, 1.48628080452,
		  1.48644006923, 1.4866038651, 1.48676986596,
		  1.48694322126, 1.48712169518, 1.48730550518,
		  1.48749488022, 1.48769006152, 1.48788829872,
		  1.48809577416, 1.48830985709, 1.48853084622,
		  1.48875905721, 1.4889948239, 1.48923485701,
		  1.4894866911, 1.48974719871, 1.49001679995,
		  1.49029594066, 1.49058509443, 1.49088027903,
		  1.49119083472, 1.49151300404, 1.49184739416,
		  1.49219465268, 1.49255547099, 1.49292496422,
		  1.49331494227, 1.49372084264, 1.49414357146,
		  1.49458410095, 1.49504347546, 1.49551561859,
		  1.49601581792, 1.4965384782, 1.49708500479,
		  1.49765691658, 1.49825585733, 1.49887415908,
		  1.49953218737, 1.5002230293, 1.50094897045,
		  1.50171250337, 1.50251635108, 1.50335070686,
		  1.50424370259, 1.50518679644};

  assert(sizeof(ope) == sizeof(ri));
  n = sizeof(ope) / sizeof(G4double);

  // Add entries into properties table
  G4MaterialPropertiesTable* pmmaMpt = new G4MaterialPropertiesTable();
  pmmaMpt -> AddProperty("RINDEX",ope, ri, n);

  fPMMA->SetMaterialPropertiesTable(pmmaMpt);

  //--------------------------------------------------
  //  Polystyrene
  //--------------------------------------------------

  G4double psOpe[] = {1.17855691476 * eV, 1.18531727948 * eV, 1.1921556484 * eV,
		  1.19907337943 * eV, 1.20724622622 * eV, 1.21434071923 * eV,
		  1.22151908801 * eV, 1.22878282887 * eV, 1.23613347391 * eV,
		  1.24407171817 * eV, 1.25173334107 * eV, 1.25961787497 * eV,
		  1.26747278096 * eV, 1.27555748388 * eV, 1.28361308037 * eV,
		  1.29190567295 * eV, 1.30016975077 * eV, 1.30867835585 * eV,
		  1.31715911434 * eV, 1.32589228353 * eV, 1.33459835773 * eV,
		  1.34356510006 * eV, 1.35250558998 * eV, 1.36171540289 * eV,
		  1.37089990528 * eV, 1.38036280821 * eV, 1.38995725822 * eV,
		  1.39952802159 * eV, 1.40939169527 * eV, 1.41923291476 * eV,
		  1.42937730497 * eV, 1.43950060877 * eV, 1.44993787198 * eV,
		  1.46035556458 * eV, 1.47109856945 * eV, 1.48182368152 * eV,
		  1.49288606181 * eV, 1.50393240457 * eV, 1.51532861688 * eV,
		  1.52671084144 * eV, 1.53845622823 * eV, 1.55018989039 * eV,
		  1.5623007489 * eV, 1.5744023801 * eV, 1.5868960378 * eV,
		  1.59938322282 * eV, 1.61227812007 * eV, 1.62516958229 * eV,
		  1.6384853632 * eV, 1.65180105826 * eV, 1.66555867051 * eV,
		  1.67954737785 * eV, 1.69354169421 * eV, 1.70800643936 * eV,
		  1.7224810702 * eV, 1.73744657277 * eV, 1.7524266775 * eV,
		  1.76791939873 * eV, 1.7834319251 * eV, 1.79948022399 * eV,
		  1.81555406989 * eV, 1.83218837643 * eV, 1.84885456954 * eV,
		  1.86610757726 * eV, 1.88339947491 * eV, 1.90130635536 * eV,
		  1.91925986739 * eV, 1.93785850943 * eV, 1.95651234706 * eV,
		  1.97584362443 * eV, 1.9952395789 * eV, 2.01534765008 * eV,
		  2.03553090516 * eV, 2.05646355006 * eV, 2.0774830334 * eV,
		  2.09929203239 * eV, 2.12156378222 * eV, 2.14394237304 * eV,
		  2.16717684728 * eV, 2.19053334687 * eV, 2.214794345 * eV,
		  2.2391942827 * eV, 2.26455136864 * eV, 2.29006626216 * eV,
		  2.31659543036 * eV, 2.34330348579 * eV, 2.37108792184 * eV,
		  2.39907483423 * eV, 2.428205786 * eV, 2.45756565774 * eV,
		  2.48814343635 * eV, 2.51897983407 * eV, 2.55111496776 * eV,
		  2.5835421428 * eV, 2.6173567117 * eV, 2.65150101439 * eV,
		  2.68713020011 * eV, 2.72313172486 * eV, 2.76072561641 * eV,
		  2.79874012264 * eV, 2.83846582951 * eV};

  G4double refractiveIndexPS[] =
  {1.57173751929, 1.57183767048, 1.57193959753,
		  1.57204334288, 1.5721667357, 1.57227457257,
		  1.57238437091, 1.57249617883, 1.57261004597,
		  1.57273383165, 1.57285410857, 1.5729787113,
		  1.57310368043, 1.57323317667, 1.5733630868,
		  1.57349773782, 1.57363285457, 1.57377293942,
		  1.57391354665, 1.57405936364, 1.5742057649,
		  1.57435763324, 1.57451015352, 1.57466841511,
		  1.57482740267, 1.57499242406, 1.57516100018,
		  1.57533042706, 1.57550636561, 1.5756832479,
		  1.57586698728, 1.57605177235, 1.57624378432,
		  1.57643695389, 1.57663774673, 1.57683982008,
		  1.57704994202, 1.5772614796, 1.5774815227,
		  1.57770313013, 1.5779337346, 1.57816606718,
		  1.57840792615, 1.57865169386, 1.57890555882,
		  1.5791615319, 1.57942821874, 1.57969723403,
		  1.57997762981, 1.58026059786, 1.58055566848,
		  1.58085852965, 1.58116438039, 1.581483537,
		  1.5818059992, 1.58214265507, 1.58248296793,
		  1.58283844142, 1.58319796321, 1.58357370069,
		  1.58395392297, 1.58435151417, 1.58475407793,
		  1.58517527362, 1.58560198783, 1.58604871996,
		  1.58650158296, 1.58697598798, 1.587457212,
		  1.58796165774, 1.58847369742, 1.58901081421,
		  1.5895563996, 1.59012911651, 1.59071129158,
		  1.59132287869, 1.59195541991, 1.5925991642,
		  1.59327624748, 1.59396589694, 1.59469188318,
		  1.59543199227, 1.59621179927, 1.59700751229,
		  1.59784670385, 1.59870385069, 1.59960874309,
		  1.60053395396, 1.60151174382, 1.60251258739,
		  1.60357150612, 1.60465665688, 1.60580615826,
		  1.60698560039, 1.60823658988, 1.60952186729,
		  1.61088698413, 1.61229150919, 1.61378547607,
		  1.61532491356, 1.6169649732};

  assert(sizeof(refractiveIndexPS) == sizeof(psOpe));
  n = sizeof(psOpe) / sizeof(G4double);

  // Add entries into properties table
  G4MaterialPropertiesTable* mptPolystyrene = new G4MaterialPropertiesTable();
  mptPolystyrene->AddProperty("RINDEX",psOpe, refractiveIndexPS, n);
  fPolystyrene->SetMaterialPropertiesTable(mptPolystyrene);

  //--------------------------------------------------
  //  ZnS(Ag)
  //--------------------------------------------------
  G4double fZnSAgOpE[] = {1.0 * eV, 2.0 * eV, 3.0 * eV, 4.0 * eV};
  G4double fZnSAgRI[] = {2.36, 2.36, 2.36, 2.36};
  G4double fZnSAgAbsLen[] = {13.2459817 * um, 13.2459817 * um,
		  13.2459817 * um, 13.2459817 * um};
  n = 4;
  G4MaterialPropertiesTable* mptZnSAg = new G4MaterialPropertiesTable();
  mptZnSAg -> AddProperty("RINDEX", fZnSAgOpE, fZnSAgRI, n);
  mptZnSAg -> AddProperty("ABSLENGTH", fZnSAgOpE, fZnSAgAbsLen, n);

  mptZnSAg -> AddConstProperty("SCINTILLATIONYIELD", 37./keV);
  mptZnSAg -> AddConstProperty("RESOLUTIONSCALE", 0.0);
  mptZnSAg -> AddConstProperty("YIELDRATIO", 1.0);


  // Emission spectrum
  G4double emiOpE[] = {2.0664 * eV, 2.25426 * eV, 2.47968 * eV, 2.5098 * eV,
		  2.5394 * eV, 2.5704040921531197 * eV, 2.5941308991576104 * eV,
		  2.6182998205783337 * eV, 2.642923329862769 * eV, 2.668014374133618 *eV,
		  2.755204386360207 * eV, 2.810308474087411 * eV, 2.8482856156291327 * eV,
		  2.8676617082524603 * eV, 2.8873032268021346 * eV, 2.92740466050772 * eV,
		  2.947876021770012 * eV, 2.968635712064167 * eV, 3.032706986425264 * eV,
		  3.099604934655233 * eV, 3.5424056396059806 * eV};
  G4double emiSpe[] = {0.045, 0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
  	  	  	  	  	   0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01875};
  n = sizeof(emiSpe) / sizeof(G4double);
  assert (sizeof(emiOpE) == sizeof(emiSpe));
  mptZnSAg -> AddProperty("FASTCOMPONENT", emiOpE, emiSpe, n);
  mptZnSAg -> AddConstProperty("FASTTIMECONSTANT", 200. * ns);
  fZnSAg -> SetMaterialPropertiesTable(mptZnSAg);
}
