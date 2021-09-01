#ifndef __AnalysisDD4hep_H__
#define __AnalysisDD4hep_H__

#include <vector>

#include "Analysis.h"

class Histos;
class SimpleTree;
class Kinematics;

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

    // perform the analysis
    void process_event();

  protected:
    int find_electron(std::vector<Clusters*> ecal_cluster, std::vector<Clusters*> hcal_cluster, double e_threshold);
    double isolation(double cone_theta, double cone_phi, std::vector<Clusters*> cluster_container, double E_threshold);

  private:
    Histos *HS;
    SimpleTree *ST;
    Kinematics *kin, *kinTrue;
    TString infileName,outfileName,outfilePrefix;
    TFile *outFile;
    Double_t eleBeamEn = 5; // GeV
    Double_t ionBeamEn = 41; // GeV
    Double_t crossingAngle = 0; // mrad
    std::map<TString,BinSet*> binSchemes;

    Weights const* weight;
    Weights const* weightJet;

    std::map<int,int> PIDtoEnum;
    std::map<int,TString> finalStateName;
    std::map<int, TString> recMethodName;

};

#endif /* __AnalysisDD4hep_H__ */
