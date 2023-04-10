# KinematicsJets

Jet methods were originally declared in the `Kinematics` class. `KinematicsJets` will inherit from the `Kinematics` class to avoid redeclaration of methods.

## Methods

### Constructor
```c
KinematicsJets(Double_t enEleBeam_, Double_t enIonBeam_, Double_t crossAng_);
```
- requires beam energies and crossing angle, set by `Analysis` (at the
  macro level)

- calls the Kinematics base class constructor - see documentation there.

---

### Calculators

#### Calculate Jet Kinematics
```c
void CalculateJetKinematics(fastjet::PseudoJet jet);
void CalculateBreitJetKinematics(fastjet::PseudoJet jet);
```
- requires Jets (from `GetJets` or `GetBreitFrameJets`)
- calculates:
  - qTjet
  - pTjet
  - zjet
  - jperp
  - ...

---

#### Find Truth Jets Matching Reco Jets
```c
void CalculateJetResolution(fastjet::PseudoJet jet);
```

- requires that GetJets and  CalculateJetKinematics has run
- for the input reconstructed jet, find the closest truth level jet
- pass the kinematics of the matched truth level jet

---

## Variables

Most of these variables are `public` to allow for easy access
- many calculator methods will modify / determine them; because of
  this, you need to understand which calculations calculate what, and
  what variables each calculation depends on
- be careful not to 'accidentally' modify them!
```c

// jet objects
int jetAlgo;
double jetRad, jetMinPt; 

std::vector<fastjet::PseudoJet> jetsRec, jetsTrue;
std::vector<fastjet::PseudoJet> breitJetsRec, breitJetsTrue;
std::map<double, int> jetConstituents;
fastjet::ClusterSequence csRec;
fastjet::ClusterSequence csTrue;

// jet variables
Double_t zjet, pTjet, qTjet, mTjet, etajet, phijet, mjet, ejet;
Double_t deltaRjet;
int matchStatusjet;
Double_t pTmtjet, mTmtjet, etamtjet, phimtjet, mmtjet, emtjet;
std::vector<double> jperp;
std::vector<double> zhad_jet;
Double_t quarkpT; // struck quark information
```
