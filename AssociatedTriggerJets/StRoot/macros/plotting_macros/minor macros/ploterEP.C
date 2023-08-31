/*this script plots eta-phi tower distribution ihist-by-ihist*/

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TTree.h"
#include "Riostream.h"
#include "TStyle.h"

void plot(Float_t R=0.3, TString trigger="MB") {

  Double_t pi=3.1415;
  TString label;
  
  TString input=Form("root/%s/R%.1lf/histos_ep_R%.1lf.root",trigger.Data(),R,R);
  TString histoname="hjet_ep";
  TString title_x = "p_{T}^{jet} (GeV/c)";
  TString title_y = "probability";
  TString figName="phi_vs_ep";  

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
//gStyle->SetNdivisions(505,"Y");
//gStyle->SetNdivisions(505,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.03,"X");
  gStyle->SetLabelSize(0.03,"Y");
  
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(1);
  
  /*TH1F* hRefmult= (TH1F*) f->Get("refmult");
  Int_t nEvents=hRefmult->GetEntries();
*/
  TFile *f =new TFile(input,"OPEN");
  TH3D* h3d= (TH3D*) f->Get(histoname);
  

  const Int_t nPlots=5; //number of plots per canvas
  const Int_t nCan=4; //number of canvases
  const Int_t ntot=nPlots*nCan;
  
  Double_t rpang1[nCan]={0, pi/4, 4*pi/8,6*pi/8};
  Double_t rpang2[nCan]={pi/8, 3*pi/8, 5*pi/8, 7*pi/8};
  Double_t ang1[nPlots]={0, pi/8, 3*pi/8, 5*pi/8 ,7*pi/8};
  Double_t ang2[nPlots]={pi/8, 3*pi/8, 5*pi/8, 7*pi/8, pi};
  
  
   TH1 *frame = new TH1I("frame", "", 1000, -10, +30);
   frame->SetXTitle(title_x); 
   frame->SetYTitle(title_y);
   frame->SetAxisRange(1E-5,1,"Y");
  
        
for(Int_t j=0; j<nCan; j++){
  TCanvas *c1= new TCanvas(Form("c1_%i",j),"",10,10,1500,900);
  c1->SetFillColor(0);
  c1->SetFrameFillColor(0);
  c1->SetGrid();
  
   gPad->SetLogy();
   TString htitle= Form("#Psi_{EP}: %.2lf#pi-%.2lf#pi",rpang1[j]/pi,rpang2[j]/pi);
   frame->SetTitle(htitle);
   frame->DrawCopy("");
   
   TString h3dname=Form("histo3d_%i",j);
   TString hname=Form("histo2d_%i",j);
   TH3D* histo3d=(TH3D*) h3d->Clone(h3dname);
   histo3d->GetZaxis()->SetRange(histo3d->GetZaxis()->FindBin(rpang1[j]),histo3d->GetZaxis()->FindBin(rpang2[j])); 
   TH2D* histo2d=(TH2D*) histo3d->Project3D("xy")->Clone(hname);
      
   TH1D *histo[nPlots];
   
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  
   for(Int_t i=0; i<nPlots; i++){
  Int_t ihist=(j*nPlots)+i;
  hname=Form("histo_%i",ihist);
  histo[i] = histo2d->ProjectionY(hname, histo2d->GetXaxis()->FindBin(ang1[i]), histo2d->GetXaxis()->FindBin(ang2[i]));
  histo[i]->Scale(1/histo[i]->Integral());
   histo[i]->SetLineColor(i+1);
   histo[i]->SetLineWidth(2);
   histo[i]->DrawCopy("esame");
   legspectra->AddEntry(histo[i], Form("#phi-#Psi: %.2lf#pi-%.2lf#pi",ang1[i]/pi,ang2[i]/pi), "lp");

   }//i-loop
   
     legspectra->DrawClone("same");

  TString output=Form("obr/%s/%s_%i.png",trigger.Data(),figName.Data(),j);
  c1->SaveAs(output);
  }//j-loop
  //f->Close();
  }
     
   
