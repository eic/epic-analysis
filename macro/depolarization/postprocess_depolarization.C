R__LOAD_LIBRARY(Sidis-eic)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_depolarization(
    // TString infile="out/depol.5x41.root"
    TString infile="out/depol.18x275.root"
) {

  PostProcessor *P = new PostProcessor(infile);

  P->Op()->Payload(
      [&P](Histos *H) {
        P->DrawSingle(H,"depolAvsQ2","COLZ",1,true);
        P->DrawSingle(H,"depolBAvsQ2","COLZ",1,true);
        P->DrawSingle(H,"depolCAvsQ2","COLZ",1,true);
        P->DrawSingle(H,"depolVAvsQ2","COLZ",1,true);
        P->DrawSingle(H,"depolWAvsQ2","COLZ",1,true);
        P->DrawSingle(H,"epsilonVsQ2","COLZ",1,true);
      }
      );

  P->Execute();
  P->Finish(); 

};
