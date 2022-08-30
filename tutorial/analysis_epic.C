R__LOAD_LIBRARY(Sidis-eic)

//
// this is currently a script to support development of AnalysisEpic;
// when AnalysisEpic is ready, this will become the tutorial script
//
// currently testing with files produced from benchmarks:
//   repo:        physics_benchmarks, from https://eicweb.phy.anl.gov/EIC/benchmarks/physics_benchmarks
//   CI stage:    finish
//   CI job:      summary
//   CI artifact: results/dis/10on100/minQ2=1/rec-dis_10x100_minQ2=1.root
// 
// test procedure:
// 1. download this artifact from a recent pipeline, and store in `datarec/epic_test/`
// 2. run this macro
//
// if you use a different artifact, edit `tutorial/s3files.epic.config`
//
//

/* EPIC simulation example
 * - note the similarity of the macro to the fast simulation
 * - you only need to swap `AnalysisDelphes` with `AnalysisEpic` to switch
 *   between fast and full simulations
 */
void analysis_epic(
    TString  infiles="tutorial/s3files.epic.config", // list of input files
    Double_t eleBeamEn=10,                           // electron beam energy [GeV]
    Double_t ionBeamEn=100,                          // ion beam energy [GeV]
    Double_t crossingAngle=-25,                      // crossing angle [mrad]
    TString  outfilePrefix="tutorial.epic"           // output filename prefix
    )
{

  // setup analysis ========================================
  // - define `AnalysisEpic` instead of `AnalysisDelphes`
  AnalysisEpic *A = new AnalysisEpic(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  A->maxEvents = 300000; // use this to limit the number of events
  A->writeSimpleTree = true;

  // set reconstruction method and final states =============================
  // - see `Analysis` constructor for methods (or other tutorials)
  A->SetReconMethod("Ele");

  A->AddFinalState("pipTrack");
  //A->AddFinalState("pimTrack");
  //A->AddFinalState("KpTrack");
  //A->AddFinalState("KmTrack");
  //A->AddFinalState("jet"); // (TODO)


  // define cuts ====================================
  // - cuts are defined the same way as bins are defined; be mindful
  //   of what bins you are defining vs. what cuts you are defining.
  //   For example, if you define a Q2 minimum cut, and you also define
  //   Q2 bins below, you may be creating more bins than you actually
  //   need, since the Q2 minimum cut actually defines another bin;
  //   in this case, your Q2 bins effectively define a Q2 minimum.
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


  // set binning scheme ====================================
  // - see other tutorials for guidance
  // - see `Analysis` constructor for available bin variables
  A->AddBinScheme("q2");
  A->AddBinScheme("x");

  // 3x2 grid of (x,Q2) bins, equal width in logarithmic scale
  A->BinScheme("q2")->BuildBins( 2, 1,    100,  true );
  A->BinScheme("x")->BuildBins(  3, 0.01, 1,    true );



  // perform the analysis ==================================
  A->Execute();

  // for reference, here is a print out of HistosDAG
  // - it lists each node, together with its inputs and outputs, which
  //   indicate the connections between the nodes
  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");

}
