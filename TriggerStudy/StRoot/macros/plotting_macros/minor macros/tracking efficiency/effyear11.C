/****************************************************************************
Author:  Stephen Horvat
Date:    26Aug2013
This macro takes the single particle efficiencies from 200GeV 2011 for pions, 
protons, and kaons and combines them using the particle to hadron ratios from
p+p spectra at 200GeV provided by Jana Bielcikova.  The weighted average of the
single particle efficiencies is computed to determine the charged hadron 
efficiency. 
*****************************************************************************/

void effyear11(Int_t set=0){
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
  //Double_t A[6],B[6],C[6];
  
  Double_t A1[] = {0.692926,0.683778,0.693047,0.702961,0.694426,0.679123};
  Double_t B1[] = {0.29562,0.293875,0.291128,0.289762,0.29542,0.29644};
  Double_t C1[] = {2.97191,2.9139,2.67972,2.19617,2.80296,3.58274};
  TF1* effpminusC = new TF1("effpminusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpminusC->SetParameters(A1[set],B1[set],C1[set]);
  
  Double_t A2[] = {0.699628,0.689117,0.699075,0.698982,0.693529,0.683874};
  Double_t B2[] = {0.288179,0.293116,0.297055,0.310106,0.289412,0.289793};
  Double_t C2[] = {3.70325,4.19614,4.10985,4.33301,3.83259,3.669};

  TF1* effpplusC = new TF1("effpplusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpplusC->SetParameters(A2[set],B2[set],C2[set]);
  
  Double_t A3[] = {0.701682,0.720091,0.713548,0.700946,0.68799,0.703425};
  Double_t B3[] = {0.293838,0.239991,0.312728,0.288453,0.315407,0.308106};
  Double_t C3[] = {1.18335,0.834732,1.14848,1.27328,1.42608,1.1738};

  TF1* effkminusC = new TF1("effkminusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effkminusC->SetParameters(A3[set],B3[set],C3[set]);
  
  Double_t A4[] = {0.696464,0.680342,0.711189,0.70634,0.69172,0.671703};
  Double_t B4[] = {0.281856,0.316867,0.291133,0.268266,0.295947,0.314304};
  Double_t C4[] = {1.22047,1.30837,1.15931,1.12009,1.28311,1.44905};
  
  TF1* effkplusC = new TF1("effkplusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effkplusC->SetParameters(A4[set],B4[set],C4[set]);
  
  Double_t A5[] = {0.697626,0.682936,0.700343,0.705278,0.701274,0.690822};
  Double_t B5[] = {0.1054,0.115793,0.141843,0.124155,0.0997088,0.103618};
  Double_t C5[] = {1.79834,2.01822,3.03629,2.08263,1.45986,1.56628};

  TF1* effpiminusC = new TF1("effpiminusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpiminusC->SetParameters(A5[set],B5[set],C5[set]);
  
  Double_t A6[] = {0.697252,0.687065,0.712922,0.704841,0.695156,0.687138};
  Double_t B6[] = {0.0885106,0.108068,0.0845974,0.105015,0.108058,0.0900939};
  Double_t C6[] = {1.41963,1.5425,1.13622,1.40576,1.63167,1.29073};

  TF1* effpiplusC = new TF1("effpiplusC","[0]*exp(-pow([1]/x,[2]))",0.1,20.);
  effpiplusC->SetParameters(A6[set],B6[set],C6[set]);
  
  Double_t p1l[7]={0,11.2496,13.3369,0.423814,0.397589,0.231914,0.191385};
  Double_t p2l[7]={0,1.48493,1.19841,3.09346,4.00881,14.526, 21.2117};
  Double_t p3l[7]={0,-11.4475,-10.0152,-13.762,-17.0636,-51.1358,-74.3381};

  Double_t p1h[7]={0,11.4916,10.7937,3.77487,2.56875,0.978768,0.533545};
  Double_t p2h[7]={0,1.22731,1.27492,1.27857,1.60754,2.31051,2.99409};
  Double_t p3h[7]={0,-10.0128,-10.1654,-10.1423,-11.3835,-13.1034,-15.161};
  
  TF1* fpplusL = new TF1("fpplusL","[0]*pow([1]+x/[2],[3])",0.1,1.2);
  fpplusL->SetParameters(p1l[5],1.,p2l[5],p3l[5]);
  
  TF1* fpminusL = new TF1("fpminusL","[0]*pow([1]+x/[2],[3])",0.1,1.2);
  fpminusL->SetParameters(p1l[6],1.,p2l[6],p3l[6]);
  
  TF1* fkplusL = new TF1("fkplusL","[0]*pow([1]+x/[2],[3])",0.1,1.2);
  fkplusL->SetParameters(p1l[3],1.,p2l[3],p3l[3]);
  
  TF1* fkminusL = new TF1("fkminusL","[0]*pow([1]+x/[2],[3])",0.1,1.2);
  fkminusL->SetParameters(p1l[4],1.,p2l[4],p3l[4]);
  
  TF1* fpiplusL = new TF1("fpiplusL","[0]*pow([1]+x/[2],[3])",0.1,1.2);
  fpiplusL->SetParameters(p1l[1],1.,p2l[1],p3l[1]);
  
  TF1* fpiminusL = new TF1("fpiminusL","[0]*pow([1]+x/[2],[3])",0.1,1.2);
  fpiminusL->SetParameters(p1l[2],1.,p2l[2],p3l[2]);
  
  TF1* fpplusH = new TF1("fpplusH","[0]*pow([1]+x/[2],[3])",1.2,20.);
  fpplusH->SetParameters(p1h[5],1.,p2h[5],p3h[5]);
  
  TF1* fpminusH = new TF1("fpminusH","[0]*pow([1]+x/[2],[3])",1.2,20.);
  fpminusH->SetParameters(p1h[6],1.,p2h[6],p3h[6]);
  
  TF1* fkplusH = new TF1("fkplusH","[0]*pow([1]+x/[2],[3])",1.2,20.);
  fkplusH->SetParameters(p1h[3],1.,p2h[3],p3h[3]);
  
  TF1* fkminusH = new TF1("fkminusH","[0]*pow([1]+x/[2],[3])",1.2,20.);
  fkminusH->SetParameters(p1h[4],1.,p2h[4],p3h[4]);
  
  TF1* fpiplusH = new TF1("fpiplusH","[0]*pow([1]+x/[2],[3])",1.2,20.);
  fpiplusH->SetParameters(p1h[1],1.,p2h[1],p3h[1]);
  
  TF1* fpiminusH = new TF1("fpiminusH","[0]*pow([1]+x/[2],[3])",1.2,20.);
  fpiminusH->SetParameters(p1h[2],1.,p2h[2],p3h[2]);
    
  
  //construct overall hadron efficiency for central at high and low pt
  TCanvas *chL = new TCanvas("chL", "hadronic efficiency 0-10%", 800, 600);
  
  TF1* effhL = new TF1("effhL","(fpiminusL*effpiminusC+fpiplusL*effpiplusC+fkminusL*effkminusC+fkplusL*effkplusC+fpminusL*effpminusC+fpplusL*effpplusC)/(fpiminusL+fpminusL+fkminusL+fpiplusL+fpplusL+fkplusL)",0.1,1.2);
  TF1* effhH = new TF1("effhH","(fpiminusH*effpiminusC+fpiplusH*effpiplusC+fkminusH*effkminusC+fkplusH*effkplusC+fpminusH*effpminusC+fpplusH*effpplusC)/(fpiminusH+fpminusH+fkminusH+fpiplusH+fpplusH+fkplusH)",1.2,20.);
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
  legeffL->AddEntry(effhL,"#epsilon_{hadron} 0-10% p_{T} < 1.2 GeV/c","L");
  legeffL->AddEntry(effhH,"#epsilon_{hadron} 0-10% p_{T} > 1.2 GeV/c","L");
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

  TFile *outfile = new TFile(Form("eff_%i.root",set+1), "UPDATE");
  outfile->cd();
  effhL->Write();
  effhH->Write(); 
  outfile->Close();
}

// below is the function that Jana provided for the particle to hadron ratios in 200GeV p+p
/*
TF1* ppratio(Int_t particle, Double_t pt) {
  //this function returns abundancy of a given particle with transverse momentum pt in "hadrons"
  // particle labeling code: 1 piplus, 2 piminus, 3 kplus, 4 kminus, 5 proton, 6 antiproton
  //for pt<=1.2 we use low pt parametrization, above the high pt one
 
  Double_t p1l[7]={0,11.2496,13.3369,0.423814,0.397589,0.231914,0.191385};
  Double_t p2l[7]={0,1.48493,1.19841,3.09346,4.00881,14.526, 21.2117};
  Double_t p3l[7]={0,-11.4475,-10.0152,-13.762,-17.0636,-51.1358,-74.3381};

  Double_t p1h[7]={0,11.4916,10.7937,3.77487,2.56875,0.978768,0.533545};
  Double_t p2h[7]={0,1.22731,1.27492,1.27857,1.60754,2.31051,2.99409};
  Double_t p3h[7]={0,-10.0128,-10.1654,-10.1423,-11.3835,-13.1034,-15.161};

  Double_t fpart=0.0;
  Double_t fhadr=0.0;
  for(int k=1;k<7;k++) {
    if(pt<=1.2) {
      fhadr+=p1l[k]*TMath::Power(1.0+pt/p2l[k],p3l[k]);
    } else {
      fhadr+=p1h[k]*TMath::Power(1.0+pt/p2h[k],p3h[k]);
    }
  }

  if(pt<=1.2) {
    //for particles with pt<=1.2 we use low pt parametrization, above the high pt parametrization
    fpart=p1l[particle]*TMath::Power(1.0+pt/p2l[particle],p3l[particle])/fhadr;
  } else {
    fpart=p1h[particle]*TMath::Power(1.0+pt/p2h[particle],p3h[particle])/fhadr;
  }
  
  // cout << pt << "   " << fpart << endl;
  return fpart; 
}
*/
