

void plotPurities(){
    // Get purity and efficiency plots from different methods output files
    // And plot together

    // TString p1 = "out/1bin/JB_dis-18x275-xm25.canvas.root";
    // TString p2 = "out/1bin/DA_dis-18x275-xm25.canvas.root";
    // TString p3 = "out/1bin/Ele_dis-18x275-xm25.canvas.root";

    TString p1 = "out/JB_dis-10x275-xm25.canvas.root";
    TString p2 = "out/DA_dis-10x275-xm25.canvas.root";
    TString p3 = "out/Ele_dis-10x275-xm25.canvas.root";

    TFile *f1 = TFile::Open(p1);
    TFile *f2 = TFile::Open(p2);
    TFile *f3 = TFile::Open(p3);

    // You should know nbins and limits in x and Q2 already
    int nx = 6;
    int nq = 4;
    double xMin = 1e-2; double xMax = 1;
    double qMin = 1; double qMax = 1000;

    TString outName = "dis-10x275-xm25";
    // TString outName = "dis-18x275-xm25";
    const int nNames = 5;
    TString histNames[nNames] = {"z_efficiency","z_purity","z_phiH_Res","z_pT_Res","z_z_Res"};//,"z_phiH_Res",
    TString labels[nNames] = {"K^{#pm} efficiency","K^{#pm} purity","#phi_{H}","p_{T}","z"};//"#phi_{H}",
    TString header = "10x275GeV (0.2 < z < 1.0)";
    // TString header = "18x275GeV (0.2 < z < 1.0)";
    double yMin = -0.1; double yMax=1.0;

    int nbinsz = 1; //Set manually...NOTE TODO

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

    // Get Rid of border frame
    gStyle->SetFrameLineWidth(0);
    mainpad->UseCurrentStyle();

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
    bool daz = false;
    bool dapt = false;
    bool daphi = false;
    bool elez = false;
    bool elept = false;
    bool elephi = false;
    bool jbz = false;
    bool jbpt = false;
    bool jbphi = false;
    bool purity = false;
    bool efficiency = false;
    TPad *lgpad = new TPad("lgpad", "lgpad", 0.84, 0.07, 0.98, 0.21);
    TLegend *lg = new TLegend(0.01,0.01,0.99,0.99);
    lg->SetHeader(header,"C");
    lg->SetTextSize(0.08);
    if (nNames>3) lg->SetNColumns(2);
    // if (nNames>4) lg->SetNColumns(3);

    for (int i=0; i<nx; i++) {
        for (int j=0; j<nq; j++) {
            THStack *hist = new THStack();
            // TH1D *h1_ = new TH1D("h1_","h1_",nNames,0,1);
            // TH1D *h2_ = new TH1D("h2_","h2_",nNames,0,1);
            // TH1D *h3_ = new TH1D("h3_","h3_",nNames,0,1);
            
            for (int k=0; k<nNames; k++) {
                TString name; name.Form("hist__"+histNames[k]+"__bin_%d_%d",i,j);
                // std::cout<<"Getting "<<name<<std::endl;//DEBUGGING
                
                TH1D *h1_; if (f1->Get(name)!=nullptr) h1_ = (TH1D*)f1->Get(name); else h1_ = nullptr;
                TH1D *h2_; if (f2->Get(name)!=nullptr) h2_ = (TH1D*)f2->Get(name); else h2_ = nullptr;
                TH1D *h3_; if (f3->Get(name)!=nullptr) h3_ = (TH1D*)f3->Get(name); else h3_ = nullptr;

                TH1D *h1 = new TH1D("h1","",nbinsz,0,1);
                TH1D *h2 = new TH1D("h2","",nbinsz,0,1);
                TH1D *h3 = new TH1D("h3","",nbinsz,0,1);
                for (int idx=1; idx<=nbinsz; idx++) {
                if (h1_!=nullptr && h2_!=nullptr && h3_!=nullptr){
                    // h1_->SetBinError(idx,0); h2_->SetBinError(idx,0); h3_->SetBinError(idx,0);
                    // std::cout<<"\t"<<f1->Get(name)<<" "<<f2->Get(name)<<" "<<f3->Get(name)<<std::endl;
                if (h1_->GetBinContent(idx)<h2_->GetBinContent(idx) && h1_->GetBinContent(idx)<h3_->GetBinContent(idx) && (histNames[k]!="z_purity" && histNames[k]!="z_efficiency")) {
                    // std::cout<<"\tbin content: "<<h1->GetBinContent(1)<<" nbins: "<<h1->GetNbinsX()<<std::endl;
                    h1->SetBinContent(idx,h1_->GetBinContent(idx));
                    if (histNames[k]=="z_z_Res") h1->SetMarkerStyle(24);
                    if (histNames[k]=="z_pT_Res") h1->SetMarkerStyle(26);
                    if (histNames[k]=="z_phiH_Res") h1->SetMarkerStyle(32);
                    // std::cout<<"\th1: marker style, color: "<<h1->GetMarkerStyle()<<" "<<h1->GetMarkerColor()<<std::endl;//DEBUGGING
                    h1->SetMarkerColor(8);
                    h1->SetMarkerSize(1);
                    h1->GetYaxis()->SetRangeUser(-0.05,1);
                    hist->Add(h1);
                    if ( /*((i==4 && j==1) || (i==5 && j==2)) &&*/ ((histNames[k]=="z_z_Res" && !jbz) || (histNames[k]=="z_pT_Res" && !jbpt) || (histNames[k]=="z_phiH_Res" && !jbphi)) ) {//TODO: Find where these actually pop up
                        
                        if (histNames[k]=="z_z_Res") jbz=true;
                        if (histNames[k]=="z_pT_Res") jbpt=true;
                        if (histNames[k]=="z_phiH_Res") jbphi=true;
                        std::cout<<"ADDING JB ENTRY "<<labels[k]<<" "<<i<<" "<<j<<std::endl;//DEBUGGING
                        lg->AddEntry(h1,"JB "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                    }
                    // continue; //IMPORTANT!
                }
                if (h2_->GetBinContent(idx)<h1_->GetBinContent(1) && h2_->GetBinContent(idx)<h3_->GetBinContent(idx) && (histNames[k]!="z_purity" && histNames[k]!="z_efficiency")) {
                    // std::cout<<"\tbin content: "<<h2->GetBinContent(1)<<" nbins: "<<h2->GetNbinsX()<<std::endl;
                    h2->SetBinContent(idx,h2_->GetBinContent(idx));
                    if (histNames[k]=="z_z_Res") h2->SetMarkerStyle(24);
                    if (histNames[k]=="z_pT_Res") h2->SetMarkerStyle(26);
                    if (histNames[k]=="z_phiH_Res") h2->SetMarkerStyle(32);
                    // std::cout<<"\th2: marker style, color: "<<h2->GetMarkerStyle()<<" "<<h2->GetMarkerColor()<<std::endl;//DEBUGGING
                    h2->SetMarkerColor(4);
                    h2->SetMarkerSize(1);
                    h2->GetYaxis()->SetRangeUser(-0.05,1);
                    hist->Add(h2);
                    if ((histNames[k]=="z_z_Res" && !daz) || (histNames[k]=="z_pT_Res" && !dapt) || (histNames[k]=="z_phiH_Res" && !daphi)){//TODO: Find where these actually pop up
                        if (histNames[k]=="z_z_Res") daz=true;
                        if (histNames[k]=="z_pT_Res") dapt=true;
                        if (histNames[k]=="z_phiH_Res") daphi=true;
                        std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                        lg->AddEntry(h2,"DA "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                    }
                    // continue;
                }
                // else {
                    // std::cout<<"\tbin content: "<<h3->GetBinContent(1)<<" nbins: "<<h3->GetNbinsX()<<std::endl;
                    h3->SetBinContent(idx,h3_->GetBinContent(idx));
                    if (histNames[k]=="z_z_Res") h3->SetMarkerStyle(24);
                    if (histNames[k]=="z_pT_Res") h3->SetMarkerStyle(26);
                    if (histNames[k]=="z_phiH_Res") h3->SetMarkerStyle(32);
                    if (histNames[k]=="z_purity") h3->SetMarkerStyle(25);
                    if (histNames[k]=="z_efficiency") h3->SetMarkerStyle(34);
                    // std::cout<<"\th3: marker style, color: "<<h3->GetMarkerStyle()<<" "<<h3->GetMarkerColor()<<std::endl;//DEBUGGING
                    h3->SetMarkerColor(2); if (histNames[k]=="z_purity") h3->SetMarkerColor(8);
                    if (histNames[k]=="z_efficiency") h3->SetMarkerColor(9);
                    h3->SetMarkerSize(1);
                    h3->GetYaxis()->SetRangeUser(-0.05,1);
                    hist->Add(h3);
                    if ((histNames[k]=="z_z_Res" && !elez) || (histNames[k]=="z_pT_Res" && !elept) || (histNames[k]=="z_phiH_Res" && !elephi) || !purity){//TODO: Find where these actually pop up
                        if (histNames[k]=="z_z_Res") elez=true;
                        if (histNames[k]=="z_pT_Res") elept=true;
                        if (histNames[k]=="z_phiH_Res") elephi=true;
                        // if (histNames[k]=="z_purity") purity=true;
                        // if (histNames[k]=="z_efficiency") efficiency=true;
                        std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                        if (histNames[k]!="z_purity" && histNames[k]!="z_efficiency") lg->AddEntry(h3,"Ele "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                        else if (!purity && histNames[k]=="z_purity") lg->AddEntry(h3,labels[k],"p"); if (histNames[k]=="z_purity") purity=true;
                        else if (!efficiency && histNames[k]=="z_efficiency") lg->AddEntry(h3,labels[k],"p"); if (histNames[k]=="z_efficiency") efficiency=true;
                    // }
                }
                }// if (h1!=nullptr && h2!=nullptr)

                if (h1_==nullptr && h2_!=nullptr && h3_!=nullptr){
                    // h2_->SetBinError(idx,0); h3_->SetBinError(idx,0);
                    // std::cout<<"\t"<<f1->Get(name)<<" "<<f2->Get(name)<<" "<<f3->Get(name)<<std::endl;
                // if (h1->GetBinContent(1)<h2->GetBinContent(1) && h1->GetBinContent(1)<h3->GetBinContent(1) && histNames[k]!="z_purity") {
                //     // std::cout<<"\tbin content: "<<h1->GetBinContent(1)<<" nbins: "<<h1->GetNbinsX()<<std::endl;
                //     if (histNames[k]=="z_z_Res") h1->SetMarkerStyle(24);
                //     if (histNames[k]=="z_pT_Res") h1->SetMarkerStyle(26);
                //     if (histNames[k]=="z_phiH_Res") h1->SetMarkerStyle(32);
                //     // std::cout<<"\th1: marker style, color: "<<h1->GetMarkerStyle()<<" "<<h1->GetMarkerColor()<<std::endl;//DEBUGGING
                //     h1->SetMarkerColor(8);
                //     h1->SetMarkerSize(1);
                //     h1->GetYaxis()->SetRangeUser(-0.05,1);
                //     hist->Add(h1);
                //     if ( ((i==4 && j==1) || (i==5 && j==2)) && ((histNames[k]=="z_z_Res" && !jbz) || (histNames[k]=="z_pT_Res" && !jbpt) || (histNames[k]=="z_phiH_Res" && !jbphi)) ) {//TODO: Find where these actually pop up
                        
                //         if (histNames[k]=="z_z_Res") jbz=true;
                //         if (histNames[k]=="z_pT_Res") jbpt=true;
                //         if (histNames[k]=="z_phiH_Res") jbphi=true;
                //         std::cout<<"ADDING JB ENTRY "<<labels[k]<<" "<<i<<" "<<j<<std::endl;//DEBUGGING
                //         lg->AddEntry(h1,"JB "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                //     }
                //     // continue; //IMPORTANT!
                // }
                if (h2_->GetBinContent(idx)<h3_->GetBinContent(idx) && (histNames[k]!="z_purity" && histNames[k]!="z_efficiency")) {
                    // std::cout<<"\tbin content: "<<h2->GetBinContent(1)<<" nbins: "<<h2->GetNbinsX()<<std::endl;
                    h2->SetBinContent(idx,h2_->GetBinContent(idx));
                    if (histNames[k]=="z_z_Res") h2->SetMarkerStyle(24);
                    if (histNames[k]=="z_pT_Res") h2->SetMarkerStyle(26);
                    if (histNames[k]=="z_phiH_Res") h2->SetMarkerStyle(32);
                    // std::cout<<"\th2: marker style, color: "<<h2->GetMarkerStyle()<<" "<<h2->GetMarkerColor()<<std::endl;//DEBUGGING
                    h2->SetMarkerColor(4);
                    h2->SetMarkerSize(1);
                    h2->GetYaxis()->SetRangeUser(-0.05,1);
                    hist->Add(h2);
                    if ((histNames[k]=="z_z_Res" && !daz) || (histNames[k]=="z_pT_Res" && !dapt) || (histNames[k]=="z_phiH_Res" && !daphi)){//TODO: Find where these actually pop up
                        if (histNames[k]=="z_z_Res") daz=true;
                        if (histNames[k]=="z_pT_Res") dapt=true;
                        if (histNames[k]=="z_phiH_Res") daphi=true;
                        std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                        lg->AddEntry(h2,"DA "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                    }
                    // continue;
                }
                // else {
                    // std::cout<<"\tbin content: "<<h3->GetBinContent(1)<<" nbins: "<<h3->GetNbinsX()<<std::endl;
                    h3->SetBinContent(idx,h3_->GetBinContent(idx));
                    if (histNames[k]=="z_z_Res") h3->SetMarkerStyle(24);
                    if (histNames[k]=="z_pT_Res") h3->SetMarkerStyle(26);
                    if (histNames[k]=="z_phiH_Res") h3->SetMarkerStyle(32);
                    if (histNames[k]=="z_purity") h3->SetMarkerStyle(25);
                    // std::cout<<"\th3: marker style, color: "<<h3->GetMarkerStyle()<<" "<<h3->GetMarkerColor()<<std::endl;//DEBUGGING
                    h3->SetMarkerColor(2); if (histNames[k]=="z_purity") h3->SetMarkerColor(8);
                    if (histNames[k]=="z_efficiency") h3->SetMarkerColor(9);
                    h3->SetMarkerSize(1);
                    h3->GetYaxis()->SetRangeUser(-0.05,1);
                    hist->Add(h3);
                    if ((histNames[k]=="z_z_Res" && !elez) || (histNames[k]=="z_pT_Res" && !elept) || (histNames[k]=="z_phiH_Res" && !elephi) || !purity){//TODO: Find where these actually pop up
                        if (histNames[k]=="z_z_Res") elez=true;
                        if (histNames[k]=="z_pT_Res") elept=true;
                        if (histNames[k]=="z_phiH_Res") elephi=true;
                        // if (histNames[k]=="z_purity") purity=true;
                        // if (histNames[k]=="z_efficiency") efficiency=true;
                        std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                        if (histNames[k]!="z_purity" && histNames[k]!="z_efficiency") lg->AddEntry(h3,"Ele "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                        else if (!purity && histNames[k]=="z_purity") lg->AddEntry(h3,labels[k],"p"); if (histNames[k]=="z_purity") purity=true;
                        else if (!efficiency && histNames[k]=="z_efficiency") lg->AddEntry(h3,labels[k],"p"); if (histNames[k]=="z_efficiency") efficiency=true;
                    // }
                }
                }// if (h2!=nullptr && h3!=nullptr)

                if (h1_!=nullptr && h2_==nullptr && h3_!=nullptr){
                    // h1_->SetBinError(1,0); h3_->SetBinError(1,0);
                    // std::cout<<"\t"<<f1->Get(name)<<" "<<f2->Get(name)<<" "<<f3->Get(name)<<std::endl;
                if (h1_->GetBinContent(idx)<h3_->GetBinContent(idx) && (histNames[k]!="z_purity" && histNames[k]!="z_efficiency")) {
                    h1->SetBinContent(idx,h1_->GetBinContent(idx));
                    // std::cout<<"\tbin content: "<<h1->GetBinContent(1)<<" nbins: "<<h1->GetNbinsX()<<std::endl;
                    if (histNames[k]=="z_z_Res") h1->SetMarkerStyle(24);
                    if (histNames[k]=="z_pT_Res") h1->SetMarkerStyle(26);
                    if (histNames[k]=="z_phiH_Res") h1->SetMarkerStyle(32);
                    // std::cout<<"\th1: marker style, color: "<<h1->GetMarkerStyle()<<" "<<h1->GetMarkerColor()<<std::endl;//DEBUGGING
                    h1->SetMarkerColor(8);
                    h1->SetMarkerSize(1);
                    h1->GetYaxis()->SetRangeUser(-0.05,1);
                    hist->Add(h1);
                    if ( /*((i==4 && j==1) || (i==5 && j==2)) &&*/ ((histNames[k]=="z_z_Res" && !jbz) || (histNames[k]=="z_pT_Res" && !jbpt) || (histNames[k]=="z_phiH_Res" && !jbphi)) ) {//TODO: Find where these actually pop up
                        
                        if (histNames[k]=="z_z_Res") jbz=true;
                        if (histNames[k]=="z_pT_Res") jbpt=true;
                        if (histNames[k]=="z_phiH_Res") jbphi=true;
                        std::cout<<"ADDING JB ENTRY "<<labels[k]<<" "<<i<<" "<<j<<std::endl;//DEBUGGING
                        lg->AddEntry(h1,"JB "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                    }
                    // continue; //IMPORTANT!
                }
                // if (h2->GetBinContent(1)<h3->GetBinContent(1) && histNames[k]!="z_purity") {
                //     // std::cout<<"\tbin content: "<<h2->GetBinContent(1)<<" nbins: "<<h2->GetNbinsX()<<std::endl;
                //     if (histNames[k]=="z_z_Res") h2->SetMarkerStyle(24);
                //     if (histNames[k]=="z_pT_Res") h2->SetMarkerStyle(26);
                //     if (histNames[k]=="z_phiH_Res") h2->SetMarkerStyle(32);
                //     // std::cout<<"\th2: marker style, color: "<<h2->GetMarkerStyle()<<" "<<h2->GetMarkerColor()<<std::endl;//DEBUGGING
                //     h2->SetMarkerColor(4);
                //     h2->SetMarkerSize(1);
                //     h2->GetYaxis()->SetRangeUser(-0.05,1);
                //     hist->Add(h2);
                //     if ((histNames[k]=="z_z_Res" && !daz) || (histNames[k]=="z_pT_Res" && !dapt) || (histNames[k]=="z_phiH_Res" && !daphi)){//TODO: Find where these actually pop up
                //         if (histNames[k]=="z_z_Res") daz=true;
                //         if (histNames[k]=="z_pT_Res") dapt=true;
                //         if (histNames[k]=="z_phiH_Res") daphi=true;
                //         std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                //         lg->AddEntry(h2,"DA "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                //     }
                //     // continue;
                // }
                // else {
                    // std::cout<<"\tbin content: "<<h3->GetBinContent(1)<<" nbins: "<<h3->GetNbinsX()<<std::endl;
                    h3->SetBinContent(idx,h3_->GetBinContent(idx));
                    if (histNames[k]=="z_z_Res") h3->SetMarkerStyle(24);
                    if (histNames[k]=="z_pT_Res") h3->SetMarkerStyle(26);
                    if (histNames[k]=="z_phiH_Res") h3->SetMarkerStyle(32);
                    if (histNames[k]=="z_purity") h3->SetMarkerStyle(25);
                    // std::cout<<"\th3: marker style, color: "<<h3->GetMarkerStyle()<<" "<<h3->GetMarkerColor()<<std::endl;//DEBUGGING
                    h3->SetMarkerColor(2); if (histNames[k]=="z_purity") h3->SetMarkerColor(8); if (histNames[k]=="z_efficiency") h3->SetMarkerColor(8);
                    h3->SetMarkerSize(1);
                    h3->GetYaxis()->SetRangeUser(-0.05,1);
                    hist->Add(h3);
                    if ((histNames[k]=="z_z_Res" && !elez) || (histNames[k]=="z_pT_Res" && !elept) || (histNames[k]=="z_phiH_Res" && !elephi) || !purity || !efficiency){//TODO: Find where these actually pop up
                        if (histNames[k]=="z_z_Res") elez=true;
                        if (histNames[k]=="z_pT_Res") elept=true;
                        if (histNames[k]=="z_phiH_Res") elephi=true;
                        
                        std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                        if (histNames[k]!="z_purity" && histNames[k]!="z_efficiency") lg->AddEntry(h3,"Ele "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                        else if (!purity && histNames[k]=="z_purity") lg->AddEntry(h3,labels[k],"p"); if (histNames[k]=="z_purity") purity=true;
                        else if (!efficiency && histNames[k]=="z_efficiency") lg->AddEntry(h3,labels[k],"p"); if (histNames[k]=="z_efficiency") efficiency=true;
                    // }
                }
                }// if (h1!=nullptr && h3!=nullptr)
                } // idx loop
    /*
                // if (h1!=nullptr && h2!=nullptr && h3==nullptr){
                //     // std::cout<<"\t"<<f1->Get(name)<<" "<<f2->Get(name)<<" "<<f3->Get(name)<<std::endl;
                // if (h1->GetBinContent(1)<h2->GetBinContent(1) && histNames[k]!="z_purity") {
                //     // std::cout<<"\tbin content: "<<h1->GetBinContent(1)<<" nbins: "<<h1->GetNbinsX()<<std::endl;
                //     if (histNames[k]=="z_z_Res") h1->SetMarkerStyle(24);
                //     if (histNames[k]=="z_pT_Res") h1->SetMarkerStyle(26);
                //     if (histNames[k]=="z_phiH_Res") h1->SetMarkerStyle(32);
                //     // std::cout<<"\th1: marker style, color: "<<h1->GetMarkerStyle()<<" "<<h1->GetMarkerColor()<<std::endl;//DEBUGGING
                //     h1->SetMarkerColor(8);
                //     h1->SetMarkerSize(1);
                //     h1->GetYaxis()->SetRangeUser(-0.05,1);
                //     hist->Add(h1);
                //     if ( ((i==4 && j==1) || (i==5 && j==2)) && ((histNames[k]=="z_z_Res" && !jbz) || (histNames[k]=="z_pT_Res" && !jbpt) || (histNames[k]=="z_phiH_Res" && !jbphi)) ) {//TODO: Find where these actually pop up
                        
                //         if (histNames[k]=="z_z_Res") jbz=true;
                //         if (histNames[k]=="z_pT_Res") jbpt=true;
                //         if (histNames[k]=="z_phiH_Res") jbphi=true;
                //         std::cout<<"ADDING JB ENTRY "<<labels[k]<<" "<<i<<" "<<j<<std::endl;//DEBUGGING
                //         lg->AddEntry(h1,"JB "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                //     }
                //     // continue; //IMPORTANT!
                // }
                // if (h2->GetBinContent(1)<h1->GetBinContent(1) && histNames[k]!="z_purity") {
                //     // std::cout<<"\tbin content: "<<h2->GetBinContent(1)<<" nbins: "<<h2->GetNbinsX()<<std::endl;
                //     if (histNames[k]=="z_z_Res") h2->SetMarkerStyle(24);
                //     if (histNames[k]=="z_pT_Res") h2->SetMarkerStyle(26);
                //     if (histNames[k]=="z_phiH_Res") h2->SetMarkerStyle(32);
                //     // std::cout<<"\th2: marker style, color: "<<h2->GetMarkerStyle()<<" "<<h2->GetMarkerColor()<<std::endl;//DEBUGGING
                //     h2->SetMarkerColor(4);
                //     h2->SetMarkerSize(1);
                //     h2->GetYaxis()->SetRangeUser(-0.05,1);
                //     hist->Add(h2);
                //     if ((histNames[k]=="z_z_Res" && !daz) || (histNames[k]=="z_pT_Res" && !dapt) || (histNames[k]=="z_phiH_Res" && !daphi)){//TODO: Find where these actually pop up
                //         if (histNames[k]=="z_z_Res") daz=true;
                //         if (histNames[k]=="z_pT_Res") dapt=true;
                //         if (histNames[k]=="z_phiH_Res") daphi=true;
                //         std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                //         lg->AddEntry(h2,"DA "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                //     }
                //     // continue;
                // }
                // // // else {
                // //     // std::cout<<"\tbin content: "<<h3->GetBinContent(1)<<" nbins: "<<h3->GetNbinsX()<<std::endl;
                // //     if (histNames[k]=="z_z_Res") h3->SetMarkerStyle(24);
                // //     if (histNames[k]=="z_pT_Res") h3->SetMarkerStyle(26);
                // //     if (histNames[k]=="z_phiH_Res") h3->SetMarkerStyle(32);
                // //     if (histNames[k]=="z_purity") h3->SetMarkerStyle(25);
                // //     // std::cout<<"\th3: marker style, color: "<<h3->GetMarkerStyle()<<" "<<h3->GetMarkerColor()<<std::endl;//DEBUGGING
                // //     h3->SetMarkerColor(2);
                // //     h3->SetMarkerSize(1);
                // //     h3->GetYaxis()->SetRangeUser(-0.05,1);
                // //     hist->Add(h3);
                // //     if ((histNames[k]=="z_z_Res" && !elez) || (histNames[k]=="z_pT_Res" && !elept) || (histNames[k]=="z_phiH_Res" && !elephi) || !purity){//TODO: Find where these actually pop up
                // //         if (histNames[k]=="z_z_Res") elez=true;
                // //         if (histNames[k]=="z_pT_Res") elept=true;
                // //         if (histNames[k]=="z_phiH_Res") elephi=true;
                // //         if (histNames[k]=="z_purity") purity=true;
                // //         std::cout<<"ADDING ENTRY "<<labels[k]<<std::endl;//DEBUGGING
                // //         if (histNames[k]!="z_purity") lg->AddEntry(h3,"Ele "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                // //         else lg->AddEntry(h3,labels[k],"p");
                // //     // }
                // // }
                // }// if (h2!=nullptr && h2!=nullptr)
    */

                // if (h1!=nullptr) {
                //     h1->SetMarkerColor(2);
                //     h1->SetMarkerSize(1);
                //     hist->Add(h1);
                //     if (i==0 && j==0){//TODO: Find where these actually pop up
                //         lg->AddEntry(h1,"JB "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                //     }
                //     continue;
                // }
                // if (h2!=nullptr) {
                //     h2->SetMarkerColor(4);
                //     h2->SetMarkerSize(1);
                //     // h1->GetXaxis()->SetNDivision
                //     hist->Add(h2);
                //     if (i==1 && j==1){//TODO: Find where these actually pop up
                //         lg->AddEntry(h2,"DA "+labels[k],"p");//NOTE: Only grabs hists that are in 0,0 bin
                //     }
                //     continue;
                // }
            }

            mainpad->cd((nvar2-j-1)*nvar1 + i + 1);
            gPad->SetLogx(intlog1);
            gPad->SetLogy(intlog2);
            gPad->SetGridy(intgrid2);
            gPad->SetGridx(intgrid1);
            // gPad->SetGrid(0,5);//NOTE: ADDED TODO: CHECK
            TString drawStr = "";
            switch(1) {//TODO: figure out how to get THStack dimension? //can't use hist->GetHistogram()->GetDimension()
                case 1:
                if (i==0 || nbinsz!=1) drawStr = "nostack p"; //NOTE: nostackb will just throw an error, don't use. /*"ex0 p nostack"*/
                else drawStr = "nostack p a";
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
                hist->GetHistogram()->GetYaxis()->SetNdivisions(10);
                if (nbinsz==1) hist->GetHistogram()->GetXaxis()->SetNdivisions(0);
                hist->GetHistogram()->GetYaxis()->SetLabelSize(0.09);
                // mainpad->SetGrid(0,0);
                // if (i!=0) {
                    // for (int idx=0; idx<5; idx++){
                        // TF1 *f5 = new TF1("f5","0.2",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                        // f5->SetLineStyle(2);
                        // f5->SetLineColor(1);
                        // f5->Draw("SAME");
                        // TF1 *f2 = new TF1("f2","0.4",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                        // f2->SetLineStyle(2);
                        // f2->SetLineColor(1);
                        // f2->Draw("SAME");
                        // TF1 *f3 = new TF1("f3","0.6",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                        // f3->SetLineStyle(2);
                        // f3->SetLineColor(1);
                        // f3->Draw();
                        // TF1 *f4 = new TF1("f4","0.8",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                        // f4->SetLineStyle(2);
                        // f4->SetLineColor(1);
                        // f4->Draw("SAME");

                    // }
                // }
                // hist->GetHistogram()->GetYaxis()->SetRangeUser(0,1);
                // hist->Draw(drawStr);
                TF1 *f1 = new TF1("f1","0.0",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                f1->SetLineColor(1);
                f1->SetLineStyle(1);
                f1->SetLineWidth(1);
                f1->Draw("SAME");

                TF1 *f2 = new TF1("f2","0.2",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                f2->SetLineColor(1);
                f2->SetLineStyle(2);
                f2->SetLineWidth(1);
                f2->Draw("SAME");

                TF1 *f3 = new TF1("f3","0.4",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                f3->SetLineColor(1);
                f3->SetLineStyle(2);
                f3->SetLineWidth(1);
                f3->Draw("SAME");

                TF1 *f4 = new TF1("f4","0.6",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                f4->SetLineColor(1);
                f4->SetLineStyle(2);
                f4->SetLineWidth(1);
                f4->Draw("SAME");

                TF1 *f5 = new TF1("f5","0.8",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                f5->SetLineColor(1);
                f5->SetLineStyle(2);
                f5->SetLineWidth(1);
                f5->Draw("SAME");

                TF1 *f6 = new TF1("f6","1.0",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                f6->SetLineColor(1);
                f6->SetLineStyle(2);
                f6->SetLineWidth(1);
                f6->Draw("SAME");

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