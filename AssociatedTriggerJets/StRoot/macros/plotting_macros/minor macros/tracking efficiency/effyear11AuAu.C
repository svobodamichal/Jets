/****************************************************************************
Author:  Stephen Horvat
Date:    30Sep2013
This macro takes the single particle efficiencies from 200GeV 2011 for pions, 
protons, and kaons and combines them using the particle to hadron ratios from
AuAu spectra at 200GeV provided by Jana Bielcikova.  The weighted average of 
the single particle efficiencies is computed to determine the charged hadron 
efficiency. 
*****************************************************************************/

void effyear11AuAu(){
  gStyle->SetOptStat(111111);
  gStyle->SetPalette(1);
  //gStyle->SetFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetOptDate(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  
  //insert efficiency parameters for 0-10% Au+Au from 200GeV
  Double_t A,B,C;
  
  A = 0.681718;
  B = 0.297466;
  C = 3.11838;
  TF1* effpminusC = new TF1("effpminusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpminusC->SetParameters(A,B,C);
  
  A = 0.686926;
  B = 0.292949;
  C = 3.81662;

  TF1* effpplusC = new TF1("effpplusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpplusC->SetParameters(A,B,C);
  
  A = 0.677432;
  B = 0.308963;
  C = 1.35116;

  TF1* effkminusC = new TF1("effkminusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effkminusC->SetParameters(A,B,C);
  
  A = 0.669133;
  B = 0.300664;
  C = 1.42316;
  
  TF1* effkplusC = new TF1("effkplusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effkplusC->SetParameters(A,B,C);
  
  A = 0.681529;
  B = 0.138215;
  C = 2.56587;

  TF1* effpiminusC = new TF1("effpiminusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpiminusC->SetParameters(A,B,C);
  
  A = 0.681312;
  B = 0.119371;
  C = 1.84488;

  TF1* effpiplusC = new TF1("effpiplusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpiplusC->SetParameters(A,B,C);
  
  /*Double_t p1l[7]={0,11.2496,13.3369,0.423814,0.397589,0.231914,0.191385};
  Double_t p2l[7]={0,1.48493,1.19841,3.09346,4.00881,14.526, 21.2117};
  Double_t p3l[7]={0,-11.4475,-10.0152,-13.762,-17.0636,-51.1358,-74.3381};

  Double_t p1h[7]={0,11.4916,10.7937,3.77487,2.56875,0.978768,0.533545};
  Double_t p2h[7]={0,1.22731,1.27492,1.27857,1.60754,2.31051,2.99409};
  Double_t p3h[7]={0,-10.0128,-10.1654,-10.1423,-11.3835,-13.1034,-15.161};*/

  Double_t p1l[7]={0,32117,30229.8,2336.8,2202.89,53.1456,44.743};
  Double_t p2l[7]={0,1.07999,1.07965,1.07435,1.07334,0.989034,0.99638};
  Double_t p3l[7]={0,0.0963853,0.0973201,0.126937,0.127424,0.286439,0.274714};

  Double_t p1h[7]={0,0,0,0,0,311.627,227.458};
  Double_t p2h[7]={0,0,0,0,0,3.77473,3.64532};
  Double_t p3h[7]={0,0,0,0,0,-17.7066,-17.2519};

  //[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))
  
  TF1* fpplusL = new TF1("fpplusL","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",0.1,2.);
  fpplusL->SetParameters(p1l[5],1.,p2l[5],p3l[5]);
  
  TF1* fpminusL = new TF1("fpminusL","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",0.1,2.);
  fpminusL->SetParameters(p1l[6],1.,p2l[6],p3l[6]);
  
  TF1* fkplusL = new TF1("fkplusL","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",0.1,2.);
  fkplusL->SetParameters(p1l[3],1.,p2l[3],p3l[3]);
  
  TF1* fkminusL = new TF1("fkminusL","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",0.1,2.);
  fkminusL->SetParameters(p1l[4],1.,p2l[4],p3l[4]);
  
  TF1* fpiplusL = new TF1("fpiplusL","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",0.1,2.);
  fpiplusL->SetParameters(p1l[1],1.,p2l[1],p3l[1]);
  
  TF1* fpiminusL = new TF1("fpiminusL","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",0.1,2.);
  fpiminusL->SetParameters(p1l[2],1.,p2l[2],p3l[2]);
  
  TF1* fpplusH = new TF1("fpplusH","[0]*pow([1]+x/[2],[3])",2.,20.);
  fpplusH->SetParameters(p1h[5],1.,p2h[5],p3h[5]);
  
  TF1* fpminusH = new TF1("fpminusH","[0]*pow([1]+x/[2],[3])",2.,20.);
  fpminusH->SetParameters(p1h[6],1.,p2h[6],p3h[6]);
  
  TF1* fkplusH = new TF1("fkplusH","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",2.,20.);
  fkplusH->SetParameters(p1l[3],1.,p2l[3],p3l[3]);
  
  TF1* fkminusH = new TF1("fkminusH","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",2.,20.);
  fkminusH->SetParameters(p1l[4],1.,p2l[4],p3l[4]);
  
  TF1* fpiplusH = new TF1("fpiplusH","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",2.,20.);
  fpiplusH->SetParameters(p1l[1],1.,p2l[1],p3l[1]);
  
  TF1* fpiminusH = new TF1("fpiminusH","[0]*x*x*pow([1]+([2]-[1])*x/[3],-[2]/([2]-[1]))",2.,20.);
  fpiminusH->SetParameters(p1l[2],1.,p2l[2],p3l[2]);
    
  
  //construct overall hadron efficiency for central at high and low pt
  TCanvas *chL = new TCanvas("chL", "hadronic efficiency 0-10%", 800, 600);
  
  TF1* effhL = new TF1("effhL","(fpiminusL*effpiminusC+fpiplusL*effpiplusC+fkminusL*effkminusC+fkplusL*effkplusC+fpminusL*effpminusC+fpplusL*effpplusC)/(fpiminusL+fpminusL+fkminusL+fpiplusL+fpplusL+fkplusL)",0.1,2.0);
  TF1* effhH = new TF1("effhH","(fpiminusH*effpiminusC+fpiplusH*effpiplusC+fkminusH*effkminusC+fkplusH*effkplusC+fpminusH*effpminusC+fpplusH*effpplusC)/(fpiminusH+fpminusH+fkminusH+fpiplusH+fpplusH+fkplusH)",2.,20.);
  TH1D* blank = new TH1D("blank","",10,0.,10.);
  blank->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  blank->GetXaxis()->CenterTitle();
  blank->GetYaxis()->SetTitle("efficiency");
  blank->GetYaxis()->CenterTitle();
  blank->SetStats(kFALSE);
  blank->DrawCopy("C");
  effhL->Draw("sameC");
  effhH->SetLineColor(kBlue);
  effhH->Draw("sameC");
  effpiplusC->SetLineColor(kRed);
  effpiplusC->Draw("sameC");
  
  TLegend* legeffL = new TLegend(0.5,0.2,0.85,0.49);
  legeffL->AddEntry(effhL,"#epsilon_{hadron} 0-10% p_{T} < 2.0 GeV/c","L");
  legeffL->AddEntry(effhH,"#epsilon_{hadron} 0-10% p_{T} > 2.0 GeV/c","L");
  legeffL->AddEntry(effpiplusC,"#epsilon_{#pi^{+}} 0-10%","L");
  legeffL->SetBorderSize(0);
  legeffL->SetFillColor(0);
  legeffL->Draw();
  chL->Print("effhadron.png");
 
  TCanvas *csingle = new TCanvas("csingle", "single particle efficiencies in 0-10%", 800, 600);
  blank->DrawCopy("C");
  effpplusC->SetLineColor(kBlue);
  effpplusC->Draw("sameC");
  effpminusC->SetLineColor(kBlue+2);
  effpminusC->SetLineStyle(2);
  effpminusC->Draw("sameC");
  effpiplusC->Draw("sameC");
  effpiminusC->SetLineColor(kRed+2);
  effpiminusC->SetLineStyle(2);
  effpiminusC->Draw("sameC");
  effkplusC->SetLineColor(kGreen+1);
  effkplusC->Draw("sameC");
  effkminusC->SetLineColor(kGreen+3);
  effkminusC->SetLineStyle(2);
  effkminusC->Draw("sameC");
  
  TLegend* lsingle = new TLegend(0.5,0.2,0.85,0.49);
  lsingle->AddEntry(effpplusC,"#epsilon_{p^{+}} 0-10%","L");
  lsingle->AddEntry(effpminusC,"#epsilon_{p^{-}}","L");
  lsingle->AddEntry(effpiplusC,"#epsilon_{#pi^{+}}","L");
  lsingle->AddEntry(effpiminusC,"#epsilon_{#pi^{-}}","L");
  lsingle->AddEntry(effkplusC,"#epsilon_{K^{+}}","L");
  lsingle->AddEntry(effkminusC,"#epsilon_{K^{-}}","L");
  lsingle->SetBorderSize(0);
  lsingle->SetFillColor(0);
  lsingle->Draw();
  csingle->Print("effsingle.png");

  TFile *outfile = new TFile("eff.root", "UPDATE");
  outfile->cd();
  effhL->Write();
  effhH->Write(); 
  outfile->Close();
}

