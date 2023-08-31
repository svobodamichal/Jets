#include <stdio.h>
#include "iomanip.h"
#include "TString.h"
#include "TLatex.h"

void plotPhenixNew2(TString version="GPC2", TString ext="pdf"){
  
	TString outputDir=Form("../../plotting_out/obr/%s/pp",version.Data());
	TString outputDir2="obr";
	
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetTextFont(63);
  gStyle->SetTextSize(15);
  gStyle->SetLabelSize(20,"xy");
  gStyle->SetLabelFont(63,"xy");
  gStyle->SetOptTitle(0);
  

  
  double sigma_pp=30.0; //NSD from the STAR PRL 2012 30.0+/-3.5 mb for charged pion plot
  double sigma_pp_pythia=27.6373; 
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1", "PHENIX charged and neutral pions");
  gPad->SetLogy();
  
  //************************************************************************
  //                           PHENIX
  //************************************************************************
  // Plot: p9033_d1x1y2 PHENIX charged pions arxiv: 1409.1907
  //these are piplus
  double p9033_d1x1y1_xval[] = { 5.39, 6.39, 7.41, 8.44, 9.71, 11.7 };
  double p9033_d1x1y1_xerrminus[] = { 0.3899999999999997, 0.3899999999999997, 0.41000000000000014, 0.4399999999999995, 0.7100000000000009, 0.6999999999999993 };
  double p9033_d1x1y1_xerrplus[] = { 0.6100000000000003, 0.6100000000000003, 0.5899999999999999, 0.5600000000000005, 1.2899999999999991, 1.3000000000000007 };
  double p9033_d1x1y1_yval[] = { 1.75E-5, 5.01E-6, 1.56E-6, 6.19E-7, 2.14E-7, 4.83E-8 };
  double p9033_d1x1y1_yerrminus[] = { 2.4515301344262524E-6, 3.6249137920783715E-7, 1.22065556157337E-7, 5.586591089385369E-8, 2.12602916254693E-8, 8.052949770115296E-9 };
  double p9033_d1x1y1_yerrplus[] = { 2.4515301344262524E-6, 3.6249137920783715E-7, 1.22065556157337E-7, 5.586591089385369E-8, 2.12602916254693E-8, 8.052949770115296E-9 };
  double p9033_d1x1y1_ystatminus[] = { 5.0E-7, 1.5E-7, 7.0E-8, 3.9E-8, 1.6E-8, 7.1E-9 };
  double p9033_d1x1y1_ystatplus[] = { 5.0E-7, 1.5E-7, 7.0E-8, 3.9E-8, 1.6E-8, 7.1E-9 };
  int p9033_d1x1y1_numpoints = 6;


  TGraphAsymmErrors *p9033_d1x1y1 = new TGraphAsymmErrors(p9033_d1x1y1_numpoints, p9033_d1x1y1_xval, p9033_d1x1y1_yval, p9033_d1x1y1_xerrminus, p9033_d1x1y1_xerrplus, p9033_d1x1y1_yerrminus, p9033_d1x1y1_yerrplus);
  p9033_d1x1y1->SetName("/HepData/9033/d1x1y1");
  p9033_d1x1y1->SetTitle("/HepData/9033/d1x1y1");
  p9033_d1x1y1->SetMarkerStyle(23);
  //p9033_d1x1y1->Draw("AP"); //charged pions

  
  //these are piminus
  double p9033_d1x1y2_xval[] = { 5.39, 6.39, 7.41, 8.44, 9.71, 11.7 };
  double p9033_d1x1y2_xerrminus[] = { 0.3899999999999997, 0.3899999999999997, 0.41000000000000014, 0.4399999999999995, 0.7100000000000009, 0.6999999999999993 };
  double p9033_d1x1y2_xerrplus[] = { 0.6100000000000003, 0.6100000000000003, 0.5899999999999999, 0.5600000000000005, 1.2899999999999991, 1.3000000000000007 };
  double p9033_d1x1y2_yval[] = { 1.49E-5, 4.3E-6, 1.283E-6, 4.94E-7, 1.57E-7, 3.57E-8 };
  double p9033_d1x1y2_yerrminus[] = { 2.039607805437114E-6, 3.1780497164141404E-7, 1.0E-7, 4.742362280551751E-8, 1.6401219466856727E-8, 6.6211781428987395E-9 };
  double p9033_d1x1y2_yerrplus[] = { 2.039607805437114E-6, 3.1780497164141404E-7, 1.0E-7, 4.742362280551751E-8, 1.6401219466856727E-8, 6.6211781428987395E-9 };
  double p9033_d1x1y2_ystatminus[] = { 4.0E-7, 1.3E-7, 6.0E-8, 3.5E-8, 1.3E-8, 6.0E-9 };
  double p9033_d1x1y2_ystatplus[] = { 4.0E-7, 1.3E-7, 6.0E-8, 3.5E-8, 1.3E-8, 6.0E-9 };
  int p9033_d1x1y2_numpoints = 6;

  
  TGraphAsymmErrors *p9033_d1x1y2 = new TGraphAsymmErrors(p9033_d1x1y2_numpoints, p9033_d1x1y2_xval, p9033_d1x1y2_yval, p9033_d1x1y2_xerrminus, p9033_d1x1y2_xerrplus, p9033_d1x1y2_yerrminus, p9033_d1x1y2_yerrplus);
  p9033_d1x1y2->SetName("/HepData/9033/d1x1y2");
  p9033_d1x1y2->SetTitle("/HepData/9033/d1x1y2");
  p9033_d1x1y2->SetMarkerStyle(23);
  // p9033_d1x1y2->Draw("P");


 
  // Plot: p1329_d1x1y1 PHENIX pi0 from Published in PRL 91,241803; HEP-EX/0304038
  double p1329_d1x1y1_xval[] = { 1.22, 1.72, 2.22, 2.73, 3.23, 3.73, 4.23, 4.73, 5.23, 
    5.74, 6.24, 6.74, 7.45, 8.46, 9.46, 10.86, 13.25 };
  double p1329_d1x1y1_xerrminus[] = { 0.21999999999999997, 0.21999999999999997, 0.2200000000000002, 0.22999999999999998, 0.22999999999999998, 0.22999999999999998, 0.23000000000000043, 0.23000000000000043, 0.23000000000000043, 
    0.2400000000000002, 0.2400000000000002, 0.2400000000000002, 0.4500000000000002, 0.46000000000000085, 0.46000000000000085, 0.8599999999999994, 1.25 };
  double p1329_d1x1y1_xerrplus[] = { 0.28, 0.28, 0.2799999999999998, 0.27, 0.27, 0.27, 0.2699999999999996, 0.2699999999999996, 0.2699999999999996, 
    0.2599999999999998, 0.2599999999999998, 0.2599999999999998, 0.5499999999999998, 0.5399999999999991, 0.5399999999999991, 1.1400000000000006, 1.75 };
  double p1329_d1x1y1_yval[] = { 0.373, 0.0605, 0.0122, 0.00331, 9.98E-4, 3.39E-4, 1.19E-4, 4.73E-5, 2.21E-5, 
    1.11E-5, 5.0E-6, 3.0E-6, 1.08E-6, 4.85E-7, 1.64E-7, 5.07E-8, 9.76E-9 };
  double p1329_d1x1y1_yerrminus[] = { 0.0278753558721678, 0.004431392698689657, 9.183286122080701E-4, 2.6644986019887495E-4, 9.243225796225038E-5, 3.5969162041949216E-5, 1.0281627546259397E-5, 4.484528939587746E-6, 2.2176113929180647E-6, 
    1.1368155919057408E-6, 5.699561386633186E-7, 3.73894370110062E-7, 1.4467566485072743E-7, 7.8300091953969E-8, 3.643200109793587E-8, 1.2767745878580136E-8, 4.3088493276047606E-9 };
  double p1329_d1x1y1_yerrplus[] = { 0.0278753558721678, 0.004431392698689657, 9.183286122080701E-4, 2.6644986019887495E-4, 9.243225796225038E-5, 3.5969162041949216E-5, 1.0281627546259397E-5, 4.484528939587746E-6, 2.2176113929180647E-6, 
    1.1368155919057408E-6, 5.699561386633186E-7, 3.73894370110062E-7, 1.4467566485072743E-7, 7.8300091953969E-8, 3.643200109793587E-8, 1.2767745878580136E-8, 4.3088493276047606E-9 };
  double p1329_d1x1y1_ystatminus[] = { 0.005968, 0.001089, 3.0500000000000004E-4, 1.1915999999999999E-4, 5.6886E-5, 2.4747000000000002E-5, 2.856E-6, 1.9865999999999997E-6, 1.105E-6, 
    4.995E-7, 3.15E-7, 2.3100000000000004E-7, 9.504000000000001E-8, 5.8200000000000005E-8, 3.1652000000000004E-8, 1.1306100000000001E-8, 4.03088E-9 };
  double p1329_d1x1y1_ystatplus[] = { 0.005968, 0.001089, 3.0500000000000004E-4, 1.1915999999999999E-4, 5.6886E-5, 2.4747000000000002E-5, 2.856E-6, 1.9865999999999997E-6, 1.105E-6, 
    4.995E-7, 3.15E-7, 2.3100000000000004E-7, 9.504000000000001E-8, 5.8200000000000005E-8, 3.1652000000000004E-8, 1.1306100000000001E-8, 4.03088E-9 };
  int p1329_d1x1y1_numpoints = 17;


  
  TGraphAsymmErrors *p1329_d1x1y1 = new TGraphAsymmErrors(p1329_d1x1y1_numpoints, p1329_d1x1y1_xval, p1329_d1x1y1_yval, p1329_d1x1y1_xerrminus, p1329_d1x1y1_xerrplus, p1329_d1x1y1_yerrminus, p1329_d1x1y1_yerrplus);
  p1329_d1x1y1->SetName("/HepData/1329/d1x1y1");
  p1329_d1x1y1->SetTitle("/HepData/1329/d1x1y1");
  p1329_d1x1y1->SetMarkerStyle(22);
  p1329_d1x1y1->SetMarkerColor(2);
  p1329_d1x1y1->SetMarkerSize(0.7);
  p1329_d1x1y1->Draw("AP"); //neutral pions

  
  // end of only PHENIX data input
  


  //************************************************************************
  //                           STAR charged pions
  //                         PRL108(2012)072302
  //************************************************************************

  //pi plus
  double star_d1x1y1_xval[] = {3.11E+00,3.36E+00,3.61E+00,3.87E+00,4.21E+00,4.72E+00,5.22E+00,5.72E+00,6.22E+00,6.73E+00,7.41E+00,8.71E+00,1.08E+01,1.31E+01};
  
  double star_d1x1y1_xerrminus[] = {3.00e+00,3.25e+00,3.50e+00,3.75e+00,4.00e+00,4.50e+00,5.00e+00,5.50e+00,6.00e+00,6.50e+00,7.00e+00,8.00e+00,1.00e+01,1.20e+01};

  double star_d1x1y1_xerrplus[] = {3.25e+00,3.50e+00,3.75e+00,4.00e+00,4.50e+00,5.00e+00,5.50e+00,6.00e+00,6.50e+00,7.00e+00,8.00e+00,1.00e+01,1.20e+01,1.50e+01};
  double star_d1x1y1_yval[] = {3.53e-05,2.01e-05,1.21e-05,7.25e-06,3.59e-06,1.51e-06,6.81e-07,3.26e-07,1.64e-07,9.02e-08,3.97e-08,9.64e-09,1.91e-09,4.34e-10};
 
  double star_d1x1y1_yerrminus[] = {3.09e-07,1.89e-07,1.26e-07,8.53e-08,3.52e-08,1.76e-08,9.38e-09,5.44e-09,3.23e-09,2.06e-09,8.60e-10,2.08e-10,6.82e-11,2.18e-11};
  double star_d1x1y1_yerrplus[] = {3.09e-07,1.89e-07,1.26e-07,8.53e-08,3.52e-08,1.76e-08,9.38e-09,5.44e-09,3.23e-09,2.06e-09,8.60e-10,2.08e-10,6.82e-11,2.18e-11};

  int star_d1x1y1_numpoints = 14;
  
   //calculate the error bars from the values given on the STAR page 3.0-3.25 with mean at 3.11 ...
  for (int i=0;i<star_d1x1y1_numpoints;i++) {
    star_d1x1y1_xerrminus[i]=star_d1x1y1_xval[i]-star_d1x1y1_xerrminus[i];
    star_d1x1y1_xerrplus[i]=-star_d1x1y1_xval[i]+star_d1x1y1_xerrplus[i];
    star_d1x1y1_yval[i]=star_d1x1y1_yval[i]*sigma_pp; //mb cross section
  }

  //pi minus

  double star_d1x1y2_xval[] = {3.11E+00,3.36E+00,3.61E+00,3.87E+00,4.21E+00,4.72E+00,5.22E+00,5.72E+00,6.22E+00,6.73E+00,7.41E+00,8.71E+00,1.08E+01,1.31E+01};
  
  double star_d1x1y2_xerrminus[] = {3.00e+00,3.25e+00,3.50e+00,3.75e+00,4.00e+00,4.50e+00,5.00e+00,5.50e+00,6.00e+00,6.50e+00,7.00e+00,8.00e+00,1.00e+01,1.20e+01};

  double star_d1x1y2_xerrplus[] = {3.25e+00,3.50e+00,3.75e+00,4.00e+00,4.50e+00,5.00e+00,5.50e+00,6.00e+00,6.50e+00,7.00e+00,8.00e+00,1.00e+01,1.20e+01,1.50e+01};
  double star_d1x1y2_yval[] = {3.24e-05,1.86e-05,1.07e-05,6.58e-06,3.19e-06,1.32e-06,5.72e-07,2.79e-07,1.38e-07,7.50e-08,3.24e-08,7.45e-09,1.53e-09,3.29e-10};

  double star_d1x1y2_yerrminus[] = {2.80e-07,1.76e-07,1.14e-07,7.84e-08,3.13e-08,1.55e-08,8.20e-09,4.73e-09,2.77e-09,1.76e-09,6.87e-10,1.81e-10,5.91e-11,1.60e-11};
  double star_d1x1y2_yerrplus[] = {2.80e-07,1.76e-07,1.14e-07,7.84e-08,3.13e-08,1.55e-08,8.20e-09,4.73e-09,2.77e-09,1.76e-09,6.87e-10,1.81e-10,5.91e-11,1.60e-11};

  int star_d1x1y2_numpoints = 14;
  
   //calculate the error bars from the values given on the STAR page 3.0-3.25 with mean at 3.11 ...
  for (int i=0;i<star_d1x1y1_numpoints;i++) {
    star_d1x1y2_xerrminus[i]=star_d1x1y2_xval[i]-star_d1x1y2_xerrminus[i];
    star_d1x1y2_xerrplus[i]=-star_d1x1y2_xval[i]+star_d1x1y2_xerrplus[i];
    star_d1x1y2_yval[i]=star_d1x1y2_yval[i]*sigma_pp; //mb cross section
  }

  
  TGraphAsymmErrors *star_d1x1y1 = new TGraphAsymmErrors(star_d1x1y1_numpoints, star_d1x1y1_xval, star_d1x1y1_yval, star_d1x1y1_xerrminus, star_d1x1y1_xerrplus, star_d1x1y1_yerrminus, star_d1x1y1_yerrplus);
  star_d1x1y1->SetName("STAR pi plus");
  star_d1x1y1->SetTitle("STAR pi plus");
  star_d1x1y1->SetMarkerStyle(20);
  star_d1x1y1->SetMarkerSize(0.7);
  star_d1x1y1->SetMarkerColor(4);
  star_d1x1y1->Draw("P"); //STAR charged pions


  TGraphAsymmErrors *star_d1x1y2 = new TGraphAsymmErrors(star_d1x1y2_numpoints, star_d1x1y2_xval, star_d1x1y2_yval, star_d1x1y2_xerrminus, star_d1x1y2_xerrplus, star_d1x1y2_yerrminus, star_d1x1y2_yerrplus);
  star_d1x1y2->SetName("STAR pi minus");
  star_d1x1y2->SetTitle("STAR pi minus");
  star_d1x1y2->SetMarkerSize(0.7);
  star_d1x1y2->SetMarkerStyle(21);
  star_d1x1y2->Draw("P"); //STAR charged pions


  TLegend *leg1 = new TLegend(0.4,0.6,0.8,0.8);
  leg1->SetTextSize(0.04);
  leg1->AddEntry(p9033_d1x1y2,"PHENIX pi-","P");
  leg1->AddEntry(p9033_d1x1y1,"PHENIX pi+","P");
  leg1->AddEntry(p1329_d1x1y1,"PHENIX pi0","P");
  leg1->AddEntry(star_d1x1y1,"STAR pi+","P");
  leg1->AddEntry(star_d1x1y2,"STAR pi-","P");
  
  leg1->Draw("same");

  c1->Print(Form("%s/phenix-star-all-pions.%s",outputDir2.Data(),ext.Data()));

  
  
  //Drawing comparison plots
  TCanvas *c2 = new TCanvas("c2", "PHENIX charged pions compared with pi0");
  
  gPad->SetLogy();
  p1329_d1x1y1->Draw("AP"); //neutral pions
  p9033_d1x1y1->Draw("P"); //charged pions plus phenix
  p9033_d1x1y2->Draw("P"); //charged pions minus phenix
  star_d1x1y1->Draw("P"); //charged pions plus star
  star_d1x1y2->Draw("P"); //charged pions minus star


  //here come PYTHIA data

  //resulting histogram

  //TH1D *hpt_pi0_ratio=new TH1D("hpt_pi0_ratio","pi0 spectrum 1/pt weighted",2000,0,100);
  // hpt_pi0_ratio->Sumw2();

  // TH1D *hpt_pi0_final=new TH1D("hpt_pi0_final","pi0 spectrum 1/pt weighted",2000,0,100);
  //  hpt_pi0_final->Sumw2();
  
  char hname[50];	
  TString fname;

  
  //double sf;
 
  int ipt1;
  
  ipt1=1000; //minbias data from PYTHIA are pthard labeled as 1000 (run with ptmin = 0,ptamx = -1 to get true MB spectrum)

  //PYTHIA MB simulation: 
  //bin size 2000 bins in 0-100 GeV -> 0.05 GeV
  //eta bin: +/-1 unit -> 2 units in total
  // 1/2pi to get invariant yield 
  //it was this
  // fname="data/pythia_jets_6428_0--1-ptconst_0.2-30.root";
  fname="data/pythia_jets_6428_new2_mb.root";
  TFile *f1 = new TFile(fname);
  cout << fname << endl;

  TCanvas *c3 = new TCanvas("c3","PYTHIA/PHENIX pi0 spectra",600,800);
  c3->SetLeftMargin(0.5);
  c3->Draw();


  TPad *p1 = new TPad("p1","p1",0.1,0.25,0.9,1.);
  //c3->Divide(1,2,0);              
  p1->SetTicks(1,1);
  p1->SetBottomMargin(0.0);
  p1->SetLeftMargin(5.5);
  p1->Draw();

  TPad *p2 = new TPad("p2","p2",0.1,0.0,0.9,0.25);
  gPad->SetLogy(0);
  p2->SetTopMargin(0.0);
  p2->SetLeftMargin(5.5);
  p2->SetTicks(1,1);
  p2->SetBottomMargin(0.4);
  p2->Draw();

  p1->cd();
  gPad->SetLogy(1);
  //neutral pions
  TH1D *hpt_pi0_pythia =  (TH1D*)f1->Get("hpt_pi0_w");
  int nevents=hevts->GetBinContent(1);
  int ncycles=nevents/50000.0;
  double sf=1./nevents*1./0.05*1./2.;  //1/neve 1/dpt 1/deta
  sf=sf/(2.0*TMath::Pi())*sigma_pp_pythia;
  hpt_pi0_pythia->Scale(sf);
  hpt_pi0_pythia->Rebin(10);
  hpt_pi0_pythia->Scale(0.1);
  hpt_pi0_pythia->GetXaxis()->SetRange(1.,30.);
  hpt_pi0_pythia->SetMarkerStyle(26);
  hpt_pi0_pythia->SetMarkerSize(0.7);
  hpt_pi0_pythia->SetMarkerColor(2);
  hpt_pi0_pythia->GetXaxis()->SetNdivisions(503,0);
  hpt_pi0_pythia->GetYaxis()->SetLabelFont(63);
  hpt_pi0_pythia->GetYaxis()->SetLabelSize(18); //was 15
  hpt_pi0_pythia->GetYaxis()->SetTitleFont(63);
  hpt_pi0_pythia->GetYaxis()->SetTitleSize(10); //was 15
  hpt_pi0_pythia->GetYaxis()->SetTitleOffset(2.2);
  hpt_pi0_pythia->GetYaxis()->SetTitle("E d^{3}#\sigma / dp^{3} (mb GeV^{-2}#it{c}^{3})");
  /*
  hpt_pi0_pythia->GetXaxis()->SetLabelFont(63);
  hpt_pi0_pythia->GetXaxis()->SetLabelSize(15);
  hpt_pi0_pythia->GetXaxis()->SetTitleFont(63);
  hpt_pi0_pythia->GetXaxis()->SetTitleSize(15);
  hpt_pi0_pythia->GetXaxis()->SetTitleOffset(1.8);
  hpt_pi0_pythia->GetYaxis()->SetNdivisions(210,0);
  hpt_pi0_pythia->GetYaxis()->SetLabelFont(63);
  hpt_pi0_pythia->GetYaxis()->SetLabelSize(15);
  hpt_pi0_pythia->GetYaxis()->SetTitleFont(63);
  hpt_pi0_pythia->GetYaxis()->SetTitleSize(15);
  hpt_pi0_pythia->GetYaxis()->SetTitleOffset(2.5);
  */
  hpt_pi0_pythia->Draw();
  p1329_d1x1y1->Draw("P");
  
  //fit PYTHIA with a function
 
  char formula[500];
  sprintf(formula,"[0]*(([1]-1)*([1]-2))/([1]*[2]*([1]*[2]+0.140*([1]-2)))*(1+(x-0.140)/([1]*[2]))**(-[1])");
  TF1 *funTsallis = new TF1("tsallis",formula);
  funTsallis->SetParameters(10.,5.7,0.120);
  funTsallis->SetParNames("norm", "n", "C");
  funTsallis->SetLineWidth(0.5);
  funTsallis->SetLineColor(2);
  hpt_pi0_pythia->Fit(funTsallis,"L0","",1.0,10.);

  //get the ratios of pi0 yields to PYTHIA and create correspoding TGraph
  double p1329_d1x1y1_yval_r[17];
  double p1329_d1x1y1_yerrminus_r[17];
  double p1329_d1x1y1_yerrplus_r[17];
  for(int ix=0;ix<17;ix++) {
    p1329_d1x1y1_yval_r[ix]=p1329_d1x1y1_yval[ix]/funTsallis(p1329_d1x1y1_xval[ix]);
    p1329_d1x1y1_yerrminus_r[ix]=p1329_d1x1y1_yerrminus[ix]/funTsallis(p1329_d1x1y1_xval[ix]);
    p1329_d1x1y1_yerrplus_r[ix]=p1329_d1x1y1_yerrplus[ix]/funTsallis(p1329_d1x1y1_xval[ix]);
    cout << p1329_d1x1y1_yval_r[ix] << endl; //these are PHENIX pi0/PYTHIA ratios 
  }

  TGraphAsymmErrors *p1329_d1x1y1_r = new TGraphAsymmErrors(p1329_d1x1y1_numpoints, p1329_d1x1y1_xval, p1329_d1x1y1_yval_r, p1329_d1x1y1_xerrminus, p1329_d1x1y1_xerrplus, p1329_d1x1y1_yerrminus_r, p1329_d1x1y1_yerrplus_r);
  p1329_d1x1y1_r->SetName("/HepData/1329/d1x1y1 ratio to PYTHIA");
  p1329_d1x1y1_r->SetTitle("/HepData/1329/d1x1y1 ratio to PYTHIA");
  p1329_d1x1y1_r->SetMarkerStyle(26);
  p1329_d1x1y1_r->SetMarkerColor(2);
  // p1329_d1x1y1_r->Draw("P"); //neutral pions/PYTHIA 

  
  TLegend *leg3 = new TLegend(0.4,0.6,0.8,0.8,"PYTHIA/Data");
  leg3->SetNColumns(2);
  leg3->SetTextSize(0.03);
  leg3->SetBorderSize(0);
  leg3->AddEntry(hpt_pi0_pythia,"#\pi^{0}","P");
  leg3->AddEntry(p1329_d1x1y1,"PHENIX PRL","P");
  leg3->Draw("same");

  
  p2->cd();
  TH2F *framePi0=new TH2F("framePi0","Pi0",100,0.,15.,100,0.,2.0);
  framePi0->GetXaxis()->SetNdivisions(503,0);
  framePi0->GetXaxis()->SetLabelFont(63);
  framePi0->GetXaxis()->SetLabelSize(15);
  framePi0->GetXaxis()->SetTitleFont(63);
  framePi0->GetXaxis()->SetTitleSize(15);
  framePi0->GetXaxis()->SetTitleOffset(4.0);
  framePi0->GetXaxis()->SetTickSize(0.1);
  framePi0->GetYaxis()->SetNdivisions(505,0);
  framePi0->GetYaxis()->SetLabelFont(63);
  framePi0->GetYaxis()->SetLabelSize(15);
  framePi0->GetYaxis()->SetTitleFont(63);
  framePi0->GetYaxis()->SetTitleSize(15);
  framePi0->GetYaxis()->SetTitleOffset(2.0);
  framePi0->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  framePi0->GetYaxis()->SetTitle("data/PYTHIA");
  framePi0->Draw();
  
  
  p1329_d1x1y1_r->Draw("P");
  
  c3->Print(Form("%s/pi0_phenix_pythia_mb_with_kevin.%s",outputDir2.Data(),ext.Data()));

  //===========================
  //positive pions
  //===========================
  
  TCanvas *c4 = new TCanvas("c4","PYTHIA/STAR,PHENIX pi+ spectra",600,800);
  c4->SetLeftMargin(0.5);
  c4->Draw();

  TPad *p1 = new TPad("p1","p1",0.1,0.25,0.9,1.);
  p1->SetTicks(1,1);
  p1->SetBottomMargin(0.0);
  p1->SetLeftMargin(5.5);
  p1->Draw();

  TPad *p2 = new TPad("p2","p2",0.1,0.0,0.9,0.25);
  gPad->SetLogy(0);
  p2->SetTopMargin(0.0);
  p2->SetLeftMargin(5.5);
  p2->SetTicks(1,1);
  p2->SetBottomMargin(0.4);
  p2->Draw();

  p1->cd();
  gPad->SetLogy(1);

  //positive pions
  TH1D *hpt_pi1_pythia =  (TH1D*)f1->Get("hpt_pi1_w");
  int nevents=hevts->GetBinContent(1);
  int ncycles=nevents/50000.0;
  double sf=1./nevents*1./0.05*1./2.;  //1/neve 1/dpt 1/deta
  sf=sf/(2.0*TMath::Pi())*sigma_pp_pythia;
  hpt_pi1_pythia->Scale(sf);
  hpt_pi1_pythia->Rebin(10);
  hpt_pi1_pythia->Scale(0.1);
  hpt_pi1_pythia->GetXaxis()->SetRange(0,30);
  hpt_pi1_pythia->SetMarkerStyle(24);
  hpt_pi1_pythia->SetMarkerSize(0.7);
  hpt_pi1_pythia->SetMarkerColor(4);
  hpt_pi1_pythia->Draw();
  p9033_d1x1y1->Draw("P");
  star_d1x1y1->Draw("P"); //charged pions plus star

  //fit positive pions PYTHIA with a function
 
  hpt_pi1_pythia->Fit(funTsallis,"L0","",1.0,10.);

  //get the ratios of STAR pi+ yields to PYTHIA and create correspoding TGraph
  double star_d1x1y1_yval_r[14];
  double star_d1x1y1_yerrminus_r[14];
  double star_d1x1y1_yerrplus_r[14];
  for(int ix=0;ix<14;ix++) {
    star_d1x1y1_yval_r[ix]=star_d1x1y1_yval[ix]/funTsallis(star_d1x1y1_xval[ix]);
    star_d1x1y1_yerrminus_r[ix]=star_d1x1y1_yerrminus[ix]/funTsallis(star_d1x1y1_xval[ix]);
    star_d1x1y1_yerrplus_r[ix]=star_d1x1y1_yerrplus[ix]/funTsallis(star_d1x1y1_xval[ix]);
    cout << star_d1x1y1_yval_r[ix] << endl; //these are STAR pi+/PYTHIA ratios 
  }

  TGraphAsymmErrors *star_d1x1y1_r = new TGraphAsymmErrors(star_d1x1y1_numpoints, star_d1x1y1_xval, star_d1x1y1_yval_r, star_d1x1y1_xerrminus, star_d1x1y1_xerrplus, star_d1x1y1_yerrminus_r, star_d1x1y1_yerrplus_r);
  star_d1x1y1_r->SetName("STAR pi+ ratio to PYTHIA");
  star_d1x1y1_r->SetTitle("STAR pi+ ratio to PYTHIA");
  star_d1x1y1_r->SetMarkerStyle(24);
  star_d1x1y1_r->SetMarkerColor(4);
 

  
  TLegend *leg4 = new TLegend(0.4,0.6,0.8,0.8);
  leg4->SetTextSize(0.04);
  leg4->AddEntry(hpt_pi1_pythia,"PYTHIA 6.28 pi+","P");
  leg4->AddEntry(p9033_d1x1y1,"PHENIX pi+","P");
  leg4->AddEntry(star_d1x1y1,"STAR pi+","P");
  leg4->Draw("same");
  

  
  p2->cd();
  TH2F *framePiPlus=new TH2F("framePiPlus","Pi0",100,0.,15.,100,0.,2.0);
  framePiPlus->GetXaxis()->SetNdivisions(503,0);
  framePiPlus->GetXaxis()->SetLabelFont(63);
  framePiPlus->GetXaxis()->SetLabelSize(15);
  framePiPlus->GetXaxis()->SetTitleFont(63);
  framePiPlus->GetXaxis()->SetTitleSize(15);
  framePiPlus->GetXaxis()->SetTitleOffset(4.0);
  framePiPlus->GetXaxis()->SetTickSize(0.1);
  framePiPlus->GetYaxis()->SetNdivisions(405,0);
  framePiPlus->GetYaxis()->SetLabelFont(63);
  framePiPlus->GetYaxis()->SetLabelSize(15);
  framePiPlus->GetYaxis()->SetTitleFont(63);
  framePiPlus->GetYaxis()->SetTitleSize(15);
  framePiPlus->GetYaxis()->SetTitleOffset(2.0);
  framePiPlus->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  framePiPlus->GetYaxis()->SetTitle("data/PYTHIA");
  framePiPlus->Draw();
  star_d1x1y1_r->Draw("P");


  
  c4->Print(Form("%s/piplus_phenix_star_pythia_mb_with_kevin.%s",outputDir2.Data(),ext.Data()));  

  // =======================================
  // negative pions
  // =======================================

  TCanvas *c5 = new TCanvas("c5","PYTHIA/STAR,PHENIX pi- spectra",600,800);
  c5->SetLeftMargin(0.5);
  c5->Draw();

  TPad *p1 = new TPad("p1","p1",0.1,0.25,0.9,1.);
  p1->SetTicks(1,1);
  p1->SetBottomMargin(0.0);
  p1->SetLeftMargin(5.5);
  p1->Draw();

  TPad *p2 = new TPad("p2","p2",0.1,0.0,0.9,0.25);
  gPad->SetLogy(0);
  p2->SetTopMargin(0.0);
  p2->SetLeftMargin(5.5);
  p2->SetTicks(1,1);
  p2->SetBottomMargin(0.4);
  p2->Draw();

  p1->cd();
  gPad->SetLogy(1);
  
  //negative pions
  TH1D *hpt_pi2_pythia =  (TH1D*)f1->Get("hpt_pi2_w");
  int nevents=hevts->GetBinContent(1);
  int ncycles=nevents/50000.0;
  double sf=1./nevents*1./0.05*1./2.;  //1/neve 1/dpt 1/deta
  sf=sf/(2.0*TMath::Pi())*sigma_pp_pythia;
  hpt_pi2_pythia->Scale(sf);
  hpt_pi2_pythia->Rebin(10);
  hpt_pi2_pythia->Scale(0.1);
  hpt_pi2_pythia->GetXaxis()->SetRange(0,30);
  hpt_pi2_pythia->SetMarkerStyle(25);
  hpt_pi2_pythia->SetMarkerSize(0.7);
  hpt_pi2_pythia->SetMarkerColor(1);
  hpt_pi2_pythia->Draw();
  p9033_d1x1y2->Draw("P");
  star_d1x1y2->Draw("P"); //charged pions minus star

  //fit negative pions PYTHIA with a function
 
  hpt_pi2_pythia->Fit(funTsallis,"L0","",1.0,10.);

  //get the ratios of STAR pi- yields to PYTHIA and create correspoding TGraph
  double star_d1x1y2_yval_r[14];
  double star_d1x1y2_yerrminus_r[14];
  double star_d1x1y2_yerrplus_r[14];
  for(int ix=0;ix<14;ix++) {
    star_d1x1y2_yval_r[ix]=star_d1x1y2_yval[ix]/funTsallis(star_d1x1y2_xval[ix]);
    star_d1x1y2_yerrminus_r[ix]=star_d1x1y2_yerrminus[ix]/funTsallis(star_d1x1y2_xval[ix]);
    star_d1x1y2_yerrplus_r[ix]=star_d1x1y2_yerrplus[ix]/funTsallis(star_d1x1y2_xval[ix]);
    cout << star_d1x1y2_yval_r[ix] << endl; //these are STAR pi+/PYTHIA ratios 
  }

  TGraphAsymmErrors *star_d1x1y2_r = new TGraphAsymmErrors(star_d1x1y2_numpoints, star_d1x1y2_xval, star_d1x1y2_yval_r, star_d1x1y2_xerrminus, star_d1x1y2_xerrplus, star_d1x1y2_yerrminus_r, star_d1x1y2_yerrplus_r);
  star_d1x1y2_r->SetName("STAR pi- ratio to PYTHIA");
  star_d1x1y2_r->SetTitle("STAR pi- ratio to PYTHIA");
  star_d1x1y2_r->SetMarkerStyle(25);
  star_d1x1y2_r->SetMarkerColor(1);
 
  
  TLegend *leg5 = new TLegend(0.4,0.6,0.8,0.8);
  leg5->SetTextSize(0.04);
  leg5->AddEntry(hpt_pi2_pythia,"PYTHIA 6.28 pi-","P");
  leg5->AddEntry(p9033_d1x1y2,"PHENIX pi-","P");
  leg5->AddEntry(star_d1x1y2,"STAR pi-","P");
  
  leg5->Draw("same");


  p2->cd();
  TH2F *framePiMinus=new TH2F("framePiMinus","Pi0",100,0.,15.,100,0.,2.0);
  framePiMinus->GetXaxis()->SetNdivisions(503,0);
  framePiMinus->GetXaxis()->SetLabelFont(63);
  framePiMinus->GetXaxis()->SetLabelSize(15);
  framePiMinus->GetXaxis()->SetTitleFont(63);
  framePiMinus->GetXaxis()->SetTitleSize(15);
  framePiMinus->GetXaxis()->SetTitleOffset(4.0);
  framePiMinus->GetXaxis()->SetTickSize(0.1);
  framePiMinus->GetYaxis()->SetNdivisions(405,0);
  framePiMinus->GetYaxis()->SetLabelFont(63);
  framePiMinus->GetYaxis()->SetLabelSize(15);
  framePiMinus->GetYaxis()->SetTitleFont(63);
  framePiMinus->GetYaxis()->SetTitleSize(15);
  framePiMinus->GetYaxis()->SetTitleOffset(2.0);
  framePiMinus->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  framePiMinus->GetYaxis()->SetTitle("data/PYTHIA");
  framePiMinus->Draw();
  star_d1x1y2_r->Draw("P");


  c5->Print(Form("%s/piminus_phenix_star_pythia_mb_with_kevin.%s",outputDir2.Data(),ext.Data()));


  //==========================
  // plot all pions with PYTHIA
  //===========================
  TCanvas *c6 = new TCanvas("c6","PYTHIA/DATA pion spectra",650,800);
  c6->SetLeftMargin(0.2);
  c6->Draw();


  TPad *p1 = new TPad("p1","p1",0.02,0.25,0.9,1.);
  p1->SetTicks(1,1);
  p1->SetBottomMargin(0.0);
  p1->SetLeftMargin(0.12);
  p1->Draw();

  TPad *p2 = new TPad("p2","p2",0.02,0.0,0.9,0.25);
  gPad->SetLogy(0);
  p2->SetTopMargin(0.04);
  p2->SetLeftMargin(0.12);
  p2->SetTicks(1,1);
  p2->SetBottomMargin(0.4);
  p2->Draw();

  p1->cd();
  gPad->SetLogy(1);
  //neutral pions
  hpt_pi0_pythia->GetXaxis()->SetTitleSize(20);
  hpt_pi0_pythia->GetYaxis()->SetTitleSize(20);
  hpt_pi0_pythia->GetYaxis()->SetTitleOffset(2.0);
  hpt_pi0_pythia->Draw();
  p1329_d1x1y1->Draw("P");
  //positive pions
  hpt_pi1_pythia->Draw("same");
  star_d1x1y1->Draw("P"); //charged pions plus star
  //negative pions
  hpt_pi2_pythia->Draw("same");
  star_d1x1y2->Draw("P"); //charged pions minus star
    
  // TLegend *leg3 = new TLegend(0.3,0.6,0.95,0.8,"p+p #\sqrt{s}=200 GeV PYTHIA/Data");
  
  TLegend *leg3 = new TLegend(0.2,0.7,0.82,0.85,"      p+p #\sqrt{s} = 200 GeV");
  leg3->SetTextSize(0.029);
  leg3->SetNColumns(3);
  //leg3->SetTextFont(63);
  leg3->SetBorderSize(0);
  leg3->AddEntry((TObject*)0,"#\pi^{0}:"," ");
  leg3->AddEntry(hpt_pi0_pythia,"PYTHIA  ","P");
  leg3->AddEntry(p1329_d1x1y1,"PHENIX PRL91 (2003) 241803","P");
  leg3->AddEntry((TObject*)0,"#\pi^{+}:"," ");
  leg3->AddEntry(hpt_pi1_pythia,"PYTHIA","P");
  //  leg3->AddEntry(star_d1x1y1,"STAR     PRL108 (2012) 072302","P");
  leg3->AddEntry(star_d1x1y1,"STAR","P");
  leg3->AddEntry((TObject*)0,"#\pi^{-}:"," ");
  leg3->AddEntry(hpt_pi2_pythia,"PYTHIA","P");
  // leg3->AddEntry(star_d1x1y2,"STAR     PRL108 (2012) 072302","P");
  leg3->AddEntry(star_d1x1y2,"STAR","P");
  leg3->Draw("same");

  //TLatex *t1 = new TLatex(8.65,6.0e-01,"X");
  TLatex *t1 = new TLatex(8.65,6.0e-01,"}"); //#scale[2.0]{\}}
  //t1->SetTextSize(0.15);
  t1->SetTextFont(63);
  t1->SetTextSize(26);
  t1->SetLineWidth(0.1);
  t1->Draw();
  TLatex *t2 = new TLatex(8.98,6.0e-01,"PRL108 (2012) 072302");
  t2->SetTextSize(16);
  t2->SetTextFont(63);
  t2->Draw();
  p2->cd();
  TH2F *framePi0_2=new TH2F("framePi0_2","Pi0",100,0.,15.,100,0.,2.0);
  framePi0_2->GetXaxis()->SetNdivisions(503,0);
  framePi0_2->GetXaxis()->SetLabelFont(63);
  framePi0_2->GetXaxis()->SetLabelSize(18); //was 15
  framePi0_2->GetXaxis()->SetTitleFont(63);
  framePi0_2->GetXaxis()->SetTitleSize(18);
  framePi0_2->GetXaxis()->SetTitleOffset(4.0);
  framePi0_2->GetXaxis()->SetTickSize(0.1);
  framePi0_2->GetYaxis()->SetNdivisions(504,0);
  framePi0_2->GetYaxis()->SetLabelFont(63);
  framePi0_2->GetYaxis()->SetLabelSize(18); //was 15
  framePi0_2->GetYaxis()->SetTitleFont(63);
  framePi0_2->GetYaxis()->SetTitleSize(18);
  framePi0_2->GetYaxis()->SetTitleOffset(2.0);
  framePi0_2->GetXaxis()->SetTitleSize(18);
  framePi0_2->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  framePi0_2->GetYaxis()->SetTitle("data/PYTHIA");
  framePi0_2->Draw();
  
  p1329_d1x1y1_r->Draw("P");
  star_d1x1y1_r->Draw("P");
  star_d1x1y2_r->Draw("P");
  TLine *line = new TLine(0,1.,15.,1.);
  line->SetLineStyle(3);
  line->SetLineColor(1);
  line->Draw("same");
  
  c6->Print(Form("%s/pions_data_pythia_6428_spin_tune_mb.%s",outputDir2.Data(),ext.Data()));
  c6->Print(Form("%s/pions_pythia_comparison.%s",outputDir.Data(),ext.Data())); //this one is for the paper
  
}  

