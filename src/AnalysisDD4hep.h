#ifndef AnalysisDD4hep_
#define AnalysisDD4hep_

#include <vector>
#include <fstream>

#include "Analysis.h"

class Histos;
class SimpleTree;
class Kinematics;
class BinSet;

class Clusters
{
  public:
    Clusters() {}
    Clusters(double E_, double x_, double y_, double z_, double theta_, double phi_) {}
    virtual ~Clusters() {}

    double E;
    double x;
    double y;
    double z;
    double theta;
    double phi;

};

class Particles
{
  public:
    int pid;
    int charge;
    int mcID;
    TLorentzVector vecPart;
};

class AnalysisDD4hep : public Analysis
{
  public:
    AnalysisDD4hep(
        TString infileName_="",
        Double_t eleBeamEn_=5,
        Double_t ionBeamEn_=41,
        Double_t crossingAngle_=0,
        TString outfilePrefix_=""
        );
    ~AnalysisDD4hep();

    void Execute() override;

    ClassDefOverride(AnalysisDD4hep,1);
};

#endif
