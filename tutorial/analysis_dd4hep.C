R__LOAD_LIBRARY(Largex)

/* full simulation (dd4hep) usage
 * - note the similarity of the macro to the fast simulation
 * - you only need to swap `AnalysisDelphes` with `AnalysisDD4hep` to switch
 *   between fast and full simulations
 * - some settings are specific to the full simulations, e.g. electron
 *   energy threshold
 * - this tutorial accesses files on S3:
 *   - alternatively, use your preferred method to run (download, mirror, etc.)
 *   - for S3, you must know the username and password, and have them in your environment:
 *     - `export S3_ACCESS_KEY=<login>`
 *     - `export S3_SECRET_KEY=<password>`
 *   - a sample list of files is in `s3files.list`, but if these files are moved on S3,
 *     then this list becomes out of date; in that case, use this list as a template
 *     (and if you want to replace it for the sake of keeping a working tutorial 
 *     example, send a pull request with the new list)
 */
void analysis_dd4hep(
    TString infiles="tutorial/s3files.list", /* FIXME: need example list, with cross sections (see other example macros) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="tutorial.dd4hep" /* output filename prefix*/)
{

  // setup analysis ========================================
  // - define `AnalysisDD4hep` instead of `AnalysisDelphes`
  AnalysisDD4hep *A = new AnalysisDD4hep(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  A->maxEvents = 300000; // use this to limit the number of events
  A->writeSimpleTree = true;

  // Set scatt. electron cuts
  A->SetEleEnergyThreshold(eleBeamEn * 0.1);  // default is 10% of beamE
  A->SetIsoCut(0.1);        // default is 10%
  A->SetIsoConeRadius(1.0); // default is 1.0


  // set reconstruction method and final states =============================
  // - see `Analysis` constructor for methods (or other tutorials)
  A->SetReconMethod("Ele");

  A->AddFinalState("pipTrack");
  //A->AddFinalState("pimTrack");
  //A->AddFinalState("KpTrack");
  //A->AddFinalState("KmTrack");
  //A->AddFinalState("jet"); // (TODO)


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
