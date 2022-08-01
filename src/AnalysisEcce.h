#ifndef AnalysisEcce_
#define AnalysisEcce_

#include <vector>
#include <fstream>

#include "Analysis.h"

class Histos;
class SimpleTree;
class Kinematics;
class BinSet;

class ClustersEE
{
  public:
    ClustersEE() {}
    ClustersEE(double E_, double x_, double y_, double z_, double theta_, double phi_) {}
    virtual ~ClustersEE() {}

    double E;
    double x;
    double y;
    double z;
    double theta;
    double phi;

};

class ParticlesEE
{
  public:
    int pid;
    int charge;
    int mcID;
    TLorentzVector vecPart;
};

class AnalysisEcce : public Analysis
{
  public:
    AnalysisEcce(
        TString infileName_="",
        Double_t eleBeamEn_=5,
        Double_t ionBeamEn_=41,
        Double_t crossingAngle_=0,
        TString outfilePrefix_=""
        );
    ~AnalysisEcce();

    void Execute() override;

    ClassDefOverride(AnalysisEcce,1);
};


static const double pimass = 0.13957061;
static const double kmass  = 0.493677;
static const double pmass = 0.938272081;
static const double emass = 0.000511;
static const double mumass = 0.105658376;


#endif
