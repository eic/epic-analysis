R__LOAD_LIBRARY(Largex)

// full simulation tutorial
// legacy version v0, before HistosDAG update
void analysis_dd4hep_v0(
    TString infiles="tutorial/s3files.list", /* single root file or file list */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="out/tutorial.dd4hep.root" /* output filename name*/)
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
  ana->AN->maxEvents = 100000; // use this to limit the number of events

  // Set scatt. electron cuts
  ana->SetEleEnergyThreshold(eleBeamEn * 0.1);  // default is 10% of beamE
  ana->SetIsoCut(0.1);        // default is 10%
  ana->SetIsoConeRadius(1.0); // default is 1.0

  // 3D binning: 2x2 grid of (x,Q2) bins, for two z bins
  ana->AN->BinScheme("q2")->BuildBins( 2, 1,    100,  true );
  ana->AN->BinScheme("x")->BuildBins(  2, 0.05, 1,    true );
  ana->AN->BinScheme("z")->BuildBin("Range", 0.2, 0.5 );
  ana->AN->BinScheme("z")->BuildBin("Range", 0.5, 0.8 );

  // perform analysis
  ana->process_event();

}
