#ifndef __AnalysisDD4hep_H__
#define __AnalysisDD4hep_H__

#include <vector>
#include <fstream>

#include "AnalysisDelphes.h"

class Histos;
class SimpleTree;
class Kinematics;
class AnalysisDelphes;
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

class AnalysisDD4hep
{
  public:
    AnalysisDD4hep(
        Double_t eleBeamEn_=5,
        Double_t ionBeamEn_=41,
        Double_t crossingAngle_=0,
        TString outfilePrefix_=""
        );
    ~AnalysisDD4hep();

    // perform the analysis
    void process_event();

    void SetEleEnergyThreshold(double e_threshold_) { fEThreshold = e_threshold_; }
    void SetIsoConeRadius(double r_ ) { fIsoR = r_; }
    void SetIsoCut(double isocut_ ) { fIsoCut = isocut_; }

    AnalysisDelphes *AN;

  protected:
    void CheckBins(BinSet *bs, std::vector<int> &v, Double_t var);
    int find_electron(std::vector<Clusters*> ecal_cluster, std::vector<Clusters*> hcal_cluster, double e_threshold);
    double isolation(double cone_theta, double cone_phi, std::vector<Clusters*> cluster_container, double E_threshold);

  private:
    const int NBINS = 50;

    Histos *HS;
    SimpleTree *ST;
    Kinematics *kin, *kinTrue;
    TString outfileName,outfilePrefix;
    TFile *outFile;
    Double_t eleBeamEn = 5; // GeV
    Double_t ionBeamEn = 41; // GeV
    Double_t crossingAngle = 0; // mrad
    std::map<TString,BinSet*> binSchemes_;

    Weights const* weight;
    Weights const* weightJet;

    std::map<int,int> PIDtoEnum_;
    std::vector<TString> infiles;

    double fEThreshold;
    double fIsoR;
    double fIsoCut;

 public:
    // Add files to TChain
    void AddFiles(TString infilename)
    {
      infiles.clear();
      std::ifstream ifstr(infilename);
      TString fname;
      while(ifstr >> fname)
	infiles.push_back(fname);
    }

    void AddFile(TString infilename)
    {
      infiles.clear();
      infiles.push_back(infilename);
    }


};

#endif /* __AnalysisDD4hep_H__ */
