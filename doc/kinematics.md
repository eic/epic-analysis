# Kinematics

## Methods

### Constructor
```c
Kinematics(Double_t enEleBeam, Double_t enIonBeam, Double_t crossAng);
```
- requires beam energies and crossing angle, set by `Analysis` (at the
  macro level)
- sets the following hard-coded settings:
  - `mainFrame`: the Lorentz frame in which calculations are performed,
    when there is ambiguity which frame to use
    - in general, we try to use Lorentz invariant formulas where
      possible
    - for frame-dependent formulas, we have these options (defined as
      `mainFrame`):
      - lab frame
      - head-on frame (default)
  - `qComponentsMethod`: how to obtain the virtual photon 4-momentum
    components, for some of the reconstruction methods (e.g., JB, DA)
    - see `GetQWNu_*` auxiliary methods below for details
    - options:
      - quadratic (default)
      - hadronic
      - electronic
- sets beam momenta and boosts
- sets other miscellaneous variables, such as default spin and
  polarization

---


### Preparation

#### For DIS Calculators
For each event, to prepare an instance of `Kinematics` for
calculations, set the following:
- `vecElectron`: scattered electron 4-momentum
- `vecEleBeam` and `vecIonBeam`: beam 4-momenta (if taking into
  account beam effects (generator level))
- Hadronic Final State (HFS)
  - for Delphes, call `GetHFS` for reconstructed or `GetTrueHFS` for
    generated; both require several objects (see code)
    - note: `GetHFS` uses `getTrackPID` for smeared PID
  - for DD4hep/Juggler, use `AddToHFS` in particle loop, then at the
    end call `SubtractElectronFromHFS` to omit the scattered electron
    - don't forget to call `ResetHFS()` beforehand
  - HFS variables `sigmah`, `Pxh`, and `Pyh`, are needed for
    reconstruction methods that require the HFS

#### For SIDIS Single-Hadron Calculators
For each hadron, set `vecHadron`, the hadron 4-momentum
- note that `CalculateDIS` should be called beforehand

#### For Jets
Currently implemented in `AnalysisDelphes` only!
- call `GetJets`, which does the work
- added flags to choose jet algorithm, jet radius, and minimum jet pT
- alternatively, call `GetBreitFrameJets` (requires Centauro)

---


### Calculators

#### Calculate DIS Kinematics
```c
Bool_t CalculateDIS(TString recmethod);
```
This is the main method for calculating the DIS kinematics
- boosts beams to head-on frame (used only as needed)
- calls `CalculateDISby*` to do the calculations according to
  `recmethod`, which specifies the reconstruction method; this
  choice is set by `Analysis` (at the macro level); see below
  - required inputs depend on choice of reconstruction method (see
    below)
  - calculated quantities:
    - Q2
    - x
    - y
    - W
    - Nu
    - virtual photon 4-momentum `vecQ`
    - 4-momentum for W (`vecW`)
- calculate additional boosts
  - C.O.M. frame of virtual photon and ion
  - ion rest frame
- calculate depolarization factors (and epsilon and gamma)


##### Reconstruction Methods
One of the following methods will be called by `CalculateDIS`:

```c
void CalculateDISbyElectron();
```
- Electron method
- requires electron 4-momentum
- calls auxiliary method `GetQWNu_electronic()` to obtain `vecQ`,
  `vecW`, W, and Nu
- calculates Q2, x, y


```c
void CalculateDISbyJB();
```
- Jacquet-Blondel method
- requires HFS variables
  - also requires electron 4-momentum, if using `GetQWNu_electronic()`
- calculates y, Q2, and x; the calculation of y depends on `mainFrame`
  since it involves the electron beam momentum
- calculates vecQ, vecW, W, and Nu, depending on `qComponentsMethod`


```c
void CalculateDISbyDA();
```
- Double Angle method
- requires HFS variables and electron 4-momentum
- calculates Q2, x, and y, which depend on `mainFrame`
- calculates vecQ, vecW, W, and Nu, depending on `qComponentsMethod`


```c
void CalculateDISbyMixed();
```
- Mixed Variables method
- requires HFS variables and electron 4-momentum
- calculates vecQ, vecW, W, and Nu, using `GetQWNu_electronic()`
- calculate Q2
- calculate y and x, which depend on `mainFrame`, since they require
  the electron beam momentum


```c
void CalculateDISbySigma();
```
- Sigma method
- requires HFS variables and electron 4-momentum
- calculates y, Q2, and x, which depend on `mainFrame`, since they
  require the scattered electron momentum
- calculates vecQ, vecW, W, and Nu, depending on `qComponentsMethod`


```c
void CalculateDISbyeSigma();
```
- eSigma method
- first calls `CalculateDISbySigma()`
- calculate `vecQ` using scattered electron, followed by Q2
- uses x from the Sigma method
- calculates y from this Q2 and x
- finally, calculates vecQ, vecW, W, and Nu, depending on `qComponentsMethod`


##### Auxiliary Methods
- These methods calculate vecQ, vecW, W, and Nu; one of these will be
  called from `CalculateDISby*` methods, depending on the choice of
  `qComponentsMethod`:

```c
void GetQWNu_electronic();
```
- requires electron momentum
- calculates 4-momenta vecQ and vecW by subtracting the scattered
  electron momenta from beam momenta
- calculates W and Nu


```c
void GetQWNu_hadronic();
```
- requires HFS
- calculates vecQ and vecW by the HFS 4-momentum and
  subtracting the ion beam (for vecQ)
- calculates W and Nu


```c
void GetQWNu_quadratic();
```
- requires HFS, along with Q2 and y
- solves quadratic equation for vecQ, vecW, W and Nu

---


#### Calculate Single-Hadron SIDIS Kinematics

```c
void CalculateHadronKinematics();
```
- requires DIS kinematics and `vecHadron`
- calculates:
  - z
  - mX (missing mass)
  - xF
  - phiH
  - phiS (requires spin-up reference, set in constructor)
  - pT
  - qT

---


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


### Boosts

Several boost methods are available. For all of these, the first
`TLorentzVector` parameter is copied to a new `TLorentzVector`; this
copy is then boosted, and passed by reference as the second
`TLorentzVector` parameter.

```c
// boost from Lab frame Lvec to photon+ion C.o.m. frame Cvec
void BoostToComFrame(TLorentzVector Lvec, TLorentzVector &Cvec);

// boost from Lab frame Lvec to Ion rest frame Ivec
void BoostToIonFrame(TLorentzVector Lvec, TLorentzVector &Ivec);

// boost from Lab frame Lvec to ion+electron Beam c.o.m. frame Bvec
void BoostToBeamComFrame(TLorentzVector Lvec, TLorentzVector &Bvec);

// tranform from Lab frame Lvec to Head-on frame Hvec
void TransformToHeadOnFrame(TLorentzVector Lvec, TLorentzVector &Hvec);

// transform from Head-on frame Hvec back to Lab frame Lvec
void TransformBackToLabFrame(TLorentzVector Hvec, TLorentzVector &Lvec);
```

```c
// tests and validation
void ValidateHeadOnFrame(); // test head-on frame boost
```

---


### Miscellaneous Methods
```c
// convert energy,mass to momentum
static Double_t EMtoP(Double_t energy, Double_t mass);

// vector projection: returns vA projected onto vB
static TVector3 Project(TVector3 vA, TVector3 vB);

// vector rejection: returns vC projected onto plane transverse to vD
static TVector3 Reject(TVector3 vC, TVector3 vD);

// calculate angle between two planes, spanned by vectors (e.g., phiH)
static Double_t PlaneAngle(TVector3 vA, TVector3 vB, TVector3 vC, TVector3 vD);

// shift angle to the range [-PI,+PI]
static Double_t AdjAngle(Double_t ang);

// particle masses
static Double_t ElectronMass() { return 0.000511; };
static Double_t ProtonMass()   { return 0.938272; };
static Double_t KaonMass()     { return 0.493677; };
static Double_t PionMass()     { return 0.139570; };
```


---


## Variables

Most of these variables are `public` to allow for easy access
- many calculator methods will modify / determine them; because of
  this, you need to understand which calculations calculate what, and
  what variables each calculation depends on
- be careful not to 'accidentally' modify them!
```c
// DIS
Double_t W, Q2, Nu, x, y, s;

// Single-hadron SIDIS
Double_t pLab, pTlab, phiLab, etaLab, z, pT, qT, mX, xF, phiH, phiS;

// HFS variables
Double_t sigmah, Pxh, Pyh;
TLorentzVector hadronSumVec;

// depolarization
Double_t gamma,epsilon;
// - factors A,B,C,V,W from [hep-ph/0611265] using notation from [1408.5721]
Double_t depolA, depolB, depolC, depolV, depolW;
// - ratios of factors, following notation of [1807.10606] eq. 2.3 (cf. eqs. 2.2a,b)
Double_t depolP1; // for A_UT*sin(phiH+phiS) (collins), A_UT*sin(3phiH-phiS) (pretzelosity)
Double_t depolP2; // for A_LL*const
Double_t depolP3; // for twist-3 A_UT
Double_t depolP4; // for A_LL*cos(phiH)

// lab-frame 4-vectors for beams, electron, hadron
TLorentzVector vecEleBeam, vecIonBeam;
TLorentzVector vecElectron, vecW, vecQ;
TLorentzVector vecHadron;

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

Constants, set in constructor (they are not actually `const`, just
generally treated as such):
```c
// nucleon transverse spin; if you set this externally,
// it must be done before calculating phiS
Int_t tSpin; // should be +1 or -1
Int_t lSpin; // should be +1 or -1

// beam polarization
Double_t polT;
Double_t polL;
Double_t polBeam;
```
