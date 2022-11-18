R__LOAD_LIBRARY(Sidis-eic)
#include <pybind11/embed.h>
namespace py = pybind11;
// ratios of histograms with y-cut enabled to those with y-cut disabled
void analysis_testml(
    TString configFile="datarec/in.config", /* delphes tree(s) */
    TString outfilePrefix="resolutions" /* output filename prefix*/
) {
  py::scoped_interpreter guard{};
  //outfilePrefix+="_DA";
  // setup analysis ========================================
  AnalysisEcce *A = new AnalysisEcce(
      configFile,
      outfilePrefix
      );

  A->maxEvents = 100; // use this to limit the number of events
  A->writeSimpleTree = true; // write SimpleTree (for one bin)
  A->SetReconMethod("ML"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  //A->AddFinalState("KpTrack"); // kaon final state


  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  //A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  //A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


  // set binning scheme ====================================
  // z ranges
  //A->AddBinScheme("z");
  //A->BinScheme("z")->BuildBin("Min",0.2); // needed?

  // y minima
  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBin("Full"); // a bin with no y-cut

  // perform the analysis ==================================
  A->Execute();

};
