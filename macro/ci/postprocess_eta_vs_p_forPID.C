// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// plot eta vs. p in one bin, for PID coverage requests
void postprocess_eta_vs_p_forPID(TString infile="out/coverage.fastsim.root") {
  PostProcessor *P = new PostProcessor(infile);
  auto draw_payload = [&P] (Histos *H) {
    auto draw_formatting = [] (auto hist, auto canv) {
      gStyle->SetPalette(kBird);
      std::vector<std::tuple<TString,TBox*,Color_t>> limits;
      limits.push_back({
          "#splitline{dRICH}{gas}",
          new TBox(
              9.5,  // p min
              1.3,  // eta min
              50.0, // p max
              3.7   // eta max
              ),
          kRed
          });
      limits.push_back({
          "#splitline{dRICH}{aerogel}",
          new TBox(
              2.0,  // p min
              1.3,  // eta min
              10.3, // p max
              3.7   // eta max
              ),
          kBlue
          });
      for(auto [name,box,color] : limits) {
        auto box_line = (TBox*) box->Clone();
        box->SetFillColorAlpha(color,0.3);
        box->Draw();
        box_line->SetLineWidth(3);
        box_line->SetFillStyle(kFEmpty);
        box_line->SetLineColor(color);
        box_line->Draw();
        auto text = new TLatex(1.1*box->GetX1(), (box->GetY1()+box->GetY2())/2, name);
        text->SetTextColor(color);
        text->Draw();
      }
    };
    P->DrawSingle(H, "etaVsP", "COLZ", 0, false, draw_formatting);
  };
  P->Op()->Payload(draw_payload);
  P->Execute();
  P->Finish(); 
}
