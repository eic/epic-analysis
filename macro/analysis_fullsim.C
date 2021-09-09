R__LOAD_LIBRARY(Largex)

#include "Analysis.h"
#include "AnalysisDD4hep.h"

void analysis_fullsim(
    TString infiles="file.list", /* single root file or file list */
    Double_t eleBeamEn=18, /* electron beam energy [GeV] */
    Double_t ionBeamEn=275, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="test.root" /* output filename prefix*/)
{

  AnalysisDD4hep *ana = new AnalysisDD4hep(eleBeamEn,
					   ionBeamEn,
					   crossingAngle,
					   outfilePrefix);

  // Using a list of multiple files
  ana->AddFiles(infiles);

  // To run over a single root file
  // ana->AddFile(infiles);

  ana->AN->writeSimpleTree = true;

  // Set scatt. electron cuts
  ana->SetEleEnergyThreshold(eleBeamEn * 0.1);  // default is 10% of beamE
  ana->SetIsoCut(0.1);        // default is 10%
  ana->SetIsoConeRadius(1.0); // default is 1.0

  ana->process_event();

}
