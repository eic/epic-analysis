#ifndef __AnalysisDD4hep_H__
#define __AnalysisDD4hep_H__

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

    // perform the analysis
    void process_event(); // TODO: rename process_event to Execute?
    void Execute() override { process_event(); };

    void SetEleEnergyThreshold(double e_threshold_) { fEThreshold = e_threshold_; }
    void SetIsoConeRadius(double r_ ) { fIsoR = r_; }
    void SetIsoCut(double isocut_ ) { fIsoCut = isocut_; }

  protected:
    int find_electron(std::vector<Clusters*> ecal_cluster, std::vector<Clusters*> hcal_cluster, double e_threshold);
    double isolation(double cone_theta, double cone_phi, std::vector<Clusters*> cluster_container, double E_threshold);

  private:

    double fEThreshold;
    double fIsoR;
    double fIsoCut;

  ClassDefOverride(AnalysisDD4hep,1);
};

#endif /* __AnalysisDD4hep_H__ */
