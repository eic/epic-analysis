R__LOAD_LIBRARY(Largex)

/* template analysis macro
 * - runs an analysis and nothing more (e.g., it does not define bins)
 * - useful if you want to use your own analysis code
 */
void analysis_template(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="tutorial.template" /* output filename prefix*/
) {

  // setup analysis ========================================
  Analysis *A = new Analysis(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );
  //A->maxEvents = 10000; // use this to limit the number of events
  //A->writeSimpleTree = true; // write SimpleTree (for one bin)

  // set reconstruction method =============================
  // - see `Analysis` constructor for methods
  A->SetReconMethod("Ele");

  // set binning scheme ====================================
  // - see `Analysis` constructor for available bin variables
  /* do nothing -> single bin histograms */

  // final states =========================================
  // - define final states; if you define none, default sets will
  //   be defined
  //A->AddFinalState("pipTrack");
  //A->AddFinalState("pimTrack");
  //A->AddFinalState("KpTrack");
  //A->AddFinalState("KmTrack");
  //A->AddFinalState("jet");

  // perform the analysis ==================================
  A->Execute();
};
