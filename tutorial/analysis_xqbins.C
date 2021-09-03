R__LOAD_LIBRARY(Largex)

/* run in a grid of (x,Q2) 2D bins
 * - various ways to make a grid are demonstrated
 * - observe how the resulting histograms differ in each (x,Q2) bin
 */
void analysis_xqbins(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="tutorial.xqbins" /* output filename prefix*/
) {

  // setup analysis ========================================
  Analysis *A = new Analysis(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  A->maxEvents = 30000; // use this to limit the number of events

  // set binning scheme ====================================

  // tutorial switch statement: change tutorial number to try out
  // the different binning implementations
  int tutorialNum = 0;
  switch(tutorialNum) {

    case 0:
      // 1D binning in Q2, equal width in linear scale
      // - arguments of BuildBins are: (numBins, min, max, log-scale bool)
      // - alternatively: BuildBins(TAxis*, log-scale bool)
      // - log-scale is false by default
      A->BinScheme("q2")->BuildBins(  3, 1, 100);
      break;

    case 1:
      // 3x3 grid of (x,Q2) bins, equal width in logarithmic scale
      A->BinScheme("q2")->BuildBins( 3, 1,    100,  true );
      A->BinScheme("x")->BuildBins(  3, 0.05, 1,    true );
      break;

    case 2:
      // alternatively: equal width in linear scale
      A->BinScheme("q2")->BuildBins( 3, 1,    100 );
      A->BinScheme("x")->BuildBins(  3, 0.05, 1   );
      break;

    case 3:
      // custom 2x2 grid (see `CutDef` class for more cut definitions)
      // - arguments of BuildBin are (cutType, a, b), where cutType is
      //   one of the types given in `../src/CutDef.cxx`, and a and b
      //   depend on which cutType
      // - various cutTypes are exemplified here:
      A->BinScheme("q2")->BuildBin("Max",10); // Q2 < 10 GeV2
      A->BinScheme("q2")->BuildBin("Min",10); // Q2 > 10 GeV2
      A->BinScheme("x")->BuildBin("Range", 0.05, 0.2 ); // 0.05 < x < 0.2
      A->BinScheme("x")->BuildBin("CenterDelta", 0.5, 0.1 ); // |x-0.5|<0.1
      break;

    case 4:
      // overlapping Q2 bins, by specifying various Q2 minima
      // - bins are arbitrary and allowed to be overlapping
      A->BinScheme("q2")->BuildBin("Min",1);
      A->BinScheme("q2")->BuildBin("Min",10);
      A->BinScheme("q2")->BuildBin("Min",50);
      break;

    case 5:
      // 3D binning: 2x2 grid of (x,Q2) bins, for two z bins
      // - you can add more dimensions, but be careful of the curse
      //   of dimensionality
      A->BinScheme("q2")->BuildBins( 2, 1,    100,  true );
      A->BinScheme("x")->BuildBins(  2, 0.05, 1,    true );
      A->BinScheme("z")->BuildBin("Range", 0.2, 0.5 );
      A->BinScheme("z")->BuildBin("Range", 0.5, 0.8 );
      break;
  };



  // perform the analysis ==================================
  A->Execute();
};
