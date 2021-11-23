

void plotPurities(){
    // Get purity and efficiency plots from different methods output files
    // And plot together

    TString p1 = "out/JB_dis-10x275-xm25.canvas.root";
    TString p2 = "out/DA_dis-10x275-xm25.canvas.root";

    TFile *f1 = TFile::Open(p1);
    TFile *f2 = TFile::Open(p2);

    // You should know nbins and limits in x and Q2 already
    int nx = 6;
    int nq = 4;
    double xMin = 1e-2; double xMax = 1;
    double qMin = 1; double qMax = 1000;

    TString outName = "dis-10x275-xm25";
    const int nNames = 3;
    TString histNames[nNames] = {"z_z_Res","z_pT_Res","z_phiH_Res"};
    TString labels[nNames] = {"z","p_{T}","#phi_{H}"};
    TString header = "10x275GeV";
    double yMin = -0.1; double yMax=1.0;

    // borrowed from postprocessor.cxx
    // default values set for nvar1==nvar2

    // used to be input args
    TString var1name="x"; int nvar1=nx; double var1low=xMin; double var1high=xMax; bool var1log=true;
    TString var2name="Q2"; int nvar2=nq; double var2low=qMin; double var2high=qMax; bool var2log=true;
    bool intlog1=false; bool intlog2=false; bool intgrid1=false; bool intgrid2=false;


    int canvx = 933;//700;
    int canvy = 800;//600;//TODO: check new numbers are better?
    double botmargin = 0.2;
    double leftmargin = 0.2;
    double xaxisy = 0.04;
    double xaxisx1 = 0.08;
    double xaxisx2 = 0.97;
    double yaxisx = 0.04;
    double yaxisy1 = 0.085;
    double yaxisy2 = 0.97;
    if(nvar1 > nvar2){
        // different canvas sizing/axis position for unequal binning
        canvx = 1100;
        canvy = 700;
        xaxisx1 = 0.075;
        xaxisx2 = 0.975;
        yaxisy1 = 0.08;
    }
    
    TString canvN = "canv_"+outName+"_all__";
    TString histN = "stack_"+outName+"_all__";
    for (int k=0; k<nNames; k++){canvN += histNames[k]+"__";histN += histNames[k]+"__";}
    TCanvas *canv = new TCanvas(canvN,canvN, canvx, canvy);
    TPad *mainpad = new TPad("mainpad", "mainpad", 0.07, 0.07, 0.98, 0.98);

    mainpad->SetFillStyle(4000);
    mainpad->Divide(nvar1,nvar2,0,0);
    mainpad->Draw();
    TLine * lDIRC = new TLine(6,-1,6,1);
    TLine * lDIRClow = new TLine(0.5,-1,0.5,1);
    TLine * lmRICH = new TLine(2,-1,2,-4);
    TLine * lDRICH = new TLine(2.5,1,2.5,4);
    lDIRC->SetLineColor(kRed);
    lDIRClow->SetLineColor(kRed);
    lmRICH->SetLineColor(kRed);
    lDRICH->SetLineColor(kRed);
    THStack* histArray[nvar1][nvar2];
    int drawpid = 0;
    //outfile->cd("/");
    // canv->Write();

    for (int i=0; i<nx; i++) {
        for (int j=0; j<nq; j++) {
            THStack *hist = new THStack();
            TLegend *lg = new TLegend(0.05,0.05,0.95,0.95);
            lg->SetHeader(header,"C");
            lg->SetTextSize(0.15);
            if (nNames>3) lg->SetNColumns(2);

            for (int k=0; k<nNames; k++) {
                TString name; name.Form("hist__"+histNames[k]+"__%d_%d",i,j);
                TH1D *h1; if (f1->Get(name)!=nullptr) h1 = (TH1D*)f1->Get(name);
                TH1D *h2; if (f1->Get(name)!=nullptr) h2 = (TH1D*)f2->Get(name);
                if (h1->GetBinContent(1)<h2->GetBinContent(1) || h2==nullptr) {
                    h1->SetMarkerColor(2);
                    hist->Add(h1);
                    if (i==0 && j==0){//TODO: Find where these actually pop up
                        lg->AddEntry(h1,"JB "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                    }
                }
                if (h2->GetBinContent(1)<h1->GetBinContent(1) || h1==nullptr) {
                    h2->SetMarkerColor(4);
                    hist->Add(h2);
                    if (i==1 && j==1){//TODO: Find where these actually pop up
                        lg->AddEntry(h2,"DA "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                    }
                }
            }

            mainpad->cd((nvar2-j-1)*nvar1 + i + 1);
            gPad->SetLogx(intlog1);
            gPad->SetLogy(intlog2);
            gPad->SetGridy(intgrid2);
            gPad->SetGridx(intgrid1);
            TString drawStr = "";
            switch(1) {//TODO: figure out how to get THStack dimension? //can't use hist->GetHistogram()->GetDimension()
                case 1:
                drawStr = "hist p nostackb"; //NOTE: nostackb will just throw an error, don't use. /*"ex0 p nostack"*/
                break;
                case 2:
                drawStr = "COLZ";
                break;
                case 3:
                drawStr = "BOX";
                break;
            };

            if( hist->GetNhists() > 0 ) {
                hist->Draw(drawStr);
                TF1 *f1 = new TF1("f1","0",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                f1->SetLineColor(1);
                f1->SetLineWidth(1);
                f1->Draw("SAME");
                if (i==0 && j==0) {
                    mainpad->cd(nvar1*nvar2);// Bottom right corner pad
                    lg->Draw();
                    mainpad->cd((nvar2-j-1)*nvar1 + i + 1);// Return to original pad
                }
                if(drawpid){
                    lDIRClow->Draw();
                    lDIRC->Draw();
                    lmRICH->Draw();
                    lDRICH->Draw();
                }
            }

        }
    }

    // Overall stuff
    canv->cd();

    TPad *newpad1 = new TPad("newpad1","full pad",0,0,1,1);
    TPad *newpad2 = new TPad("newpad2","full pad",0,0,1,1);
    newpad1->SetFillStyle(4000);
    newpad1->Draw();
    newpad2->SetFillStyle(4000);
    newpad2->Draw();

    TString xopt, yopt;
    if(var1log) xopt = "GS";
    else xopt = "S";
    if(var2log) yopt = "GS";
    else yopt = "S";

    TGaxis *xaxis = new TGaxis(xaxisx1,xaxisy,xaxisx2,xaxisy,var1low,var1high,510,xopt);
    TGaxis *yaxis = new TGaxis(yaxisx,yaxisy1,yaxisx,yaxisy2,var2low,var2high,510,yopt);
    xaxis->SetTitle(var1name);
    xaxis->SetName("xaxis");
    xaxis->SetTitleSize(0.02);
    xaxis->SetTextFont(40);
    xaxis->SetLabelSize(0.02);
    xaxis->SetTickSize(0.02);

    yaxis->SetTitle(var2name);
    yaxis->SetTitleSize(0.02);
    yaxis->SetName("yaxis");
    yaxis->SetTextFont(40);
    yaxis->SetLabelSize(0.02);
    yaxis->SetTickSize(0.02);

    newpad1->cd();
    yaxis->Draw();
    newpad2->cd();
    xaxis->Draw();

    //  canv->Write();
    canv->Print(canvN+".png");
    canv->Print(canvN+".pdf");

    // OLD
    // // THStack *hs = new THStack("myStack","myStack");
    // // hs->Add(h1); hs->Add(h2);

    // // Make legend
    // TString header = "testheader";
    // TLegend *lg = new TLegend(0.80,0.05,0.90,0.15);
    // lg->SetHeader(header,"C");
    // lg->SetTextSize(0.15);
    // lg->SetNColumns(2);
    // lg->AddEntry(h1,"label1","p");
    // lg->AddEntry(h2,"label2","p");

    // // Draw and save
    // TCanvas *c1 = new TCanvas();
    // h1->Draw("BOX");
    // h2->Draw("SAME")
    // lg->Draw();
    // c1->Print("myStack.pdf");

}