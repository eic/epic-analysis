#include "Modulation.h"

ClassImp(Modulation)

// primary constructor
Modulation::Modulation(Int_t ID_, Int_t polarization_)
  : ID(ID_)
  , polarization(polarization_)
{ Initialize(); };


// alternative constructor, for parsing BruAsymmetry 
// formatted amplitude names
Modulation::Modulation(TString ampStr) {
  sscanf(ampStr,"AmpP%dI%d",&polarization,&ID);
  Initialize();
};


// build formula (called by the constructor)
// see [hep-ph/0611265] and [1807.10606 t.o.c.]
void Modulation::Initialize() {
  baseStr = "0";

  /* define modulation properties:
   * - `baseStr`: modulation formula string
   * - `twist`: twist of associated structure function
   * - `polT`: polarization title (which may include more than just the polarizations)
   * - `depol`: associated depolarization factor ratio p0,p1,p2,p3,p4, from [1807.10606 eq. 2.3], where p0=1
   */
  if(polarization==kUT) {
    switch(ID) {
      case  0:  baseStr="sin(phiH-phiS)";    twist=2;  polT="UT,T";  depol=0;  break;
      case  1:  baseStr="sin(phiH+phiS)";    twist=2;  polT="UT";    depol=1;  break;
      case  2:  baseStr="sin(3*phiH-phiS)";  twist=2;  polT="UT";    depol=1;  break;
      case  3:  baseStr="sin(phiS)";         twist=3;  polT="UT";    depol=3;  break;
      case  4:  baseStr="sin(2*phiH-phiS)";  twist=3;  polT="UT";    depol=3;  break;
    };
  }
  //else if(polarization==kLL) { ... };
  if(baseStr=="0") fprintf(stderr,"ERROR: bad modulation ID or polarization\n");
};


// build formula string for BruFit
TString Modulation::FormuBru() {
  formuStr = baseStr;
  Tools::GlobalRegexp(formuStr,TRegexp("phiH"),"@PhiH[]");
  Tools::GlobalRegexp(formuStr,TRegexp("phiS"),"@PhiS[]");
  return formuStr;
};


// amplitude name for parameter name
TString Modulation::AmpName() {
  TString retstr = Form("AmpP%dI%d",polarization,ID);
  return retstr;
};


// modulation formula for ROOT title
TString Modulation::ModulationTitle() {
  TString retstr = baseStr;
  retstr.ReplaceAll("phi","#phi");
  return retstr;
};


// asymmetry title
TString Modulation::AsymmetryTitle() {
  TString retstr = 
    "A_{"+PolarizationTitle()+"}["+ ModulationTitle()+"]";
  return retstr;
};


Modulation::~Modulation() {};
