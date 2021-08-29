R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

void postprocess_pTvsEta(
     TString infile = "out/coverage_20x10.cross_5x41_25_newcard.root"
){
  Analysis *A = new Analysis();
  PostProcessor *P = new PostProcessor(infile);

  // number of bins in x and Q2
  int nx = 20;
  int nq2 = 10;

  // just counters for filling Histos name vector
  int xbin = 0;
  int q2bin = 0;
  
  // initialize this 2D vector to be some large size
  std::vector<std::vector<TString>> histNames_xQ2(30,std::vector<TString>(30));

  
  for(int by  : P->GetBinNums("y")) {
  if(P->GetBinCut("y",by)->GetCutType()=="Full"){
  for(int bpt : P->GetBinNums("pt")) {
  if(P->GetBinCut("pt",bpt)->GetCutType()=="Full"){    
  for(int bz  : P->GetBinNums("z")) {
  if(P->GetBinCut("z",bz)->GetCutType()=="Full"){
  for(int bfs : P->GetBinNums("finalState")) {
    if(P->GetBinCut("finalState",bfs)->GetCutTitle() =="#pi^{+} tracks"){ // is there a better way to select a type of track?
    for(int bx  : P->GetBinNums("x")) {
      if(P->GetBinCut("x",bx)->GetCutType()=="Range"){  // not perfect, as there could be some case with another range outside of those we want	
	for(int bq  : P->GetBinNums("q2")) {
	  if(P->GetBinCut("q2",bq)->GetCutType()=="Range"){	
	    histNames_xQ2[xbin][q2bin] = A->GetHistosName(bpt,bx,bz,bq,by,bfs);	
	    q2bin++;
	  };
	};
	xbin++;
	q2bin=0;
      }; 
    };

    P->DrawInBins("xQ2_piplus", histNames_xQ2, "etaVsP", "x", nx, 1e-3, 1, true, "Q^{2}", nq2, 9.99, 1000, true);
    xbin = 0;
    q2bin = 0;

    }}}}}}}};
  
  P->Finish();
};
