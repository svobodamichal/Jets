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

void plot(Float_t R=0.3, Int_t pTlead=5, TString trigger="MB") {

  Double_t pi=3.1415;
  TString label;
  
  TString input=Form("root/%s/R%.1lf/histos_ep_R%.1lf.root",trigger.Data(),R,R);
  TString histoname=Form("hjet_ep%i",pTlead);
  TString title_x = "p_{T}^{jet} (GeV/c)";
  TString title_y = "probability";
  TString figName="sumpT_phi_vs_ep";  

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
  TFile *fout =new TFile(Form("histo_measured_pTlead%i.root",pTlead),"RECREATE");
  TFile *f =new TFile(input,"OPEN");
  TH3D* h3d= (TH3D*) f->Get(histoname);
  TH2D* histold2d=(TH2D*) f->Get("hpT_pTlead_nobadsecedge1R");
  TH1D* histoold1=histold2d->ProjectionX("histoold1",histold2d->GetYaxis()->FindBin(pTlead),histold2d->GetYaxis()->GetNbins());
  TH1D* histoold2=(TH1D*) h3d->ProjectionX("histoold2",0,-1,0,-1);
  TH1D* histosum;
   
  const Int_t nPlots=5; //number of plots per canvas
  const Int_t nCan=4; //number of canvases
  const Int_t ntot=nPlots*nCan;
  
  Double_t rpang1[nCan]={0, pi/4, pi/2,3*pi/4};
  Double_t rpang2[nCan]={pi/4-0.001, pi/2, 3*pi/4, pi};
  Double_t ang1[nPlots]={0, pi/8, 3*pi/8, 5*pi/8 ,7*pi/8};
  Double_t ang2[nPlots]={pi/8, 3*pi/8, 5*pi/8, 7*pi/8, pi};
  
  
   TH1I *frame = new TH1I("frame", "", 1000, -10, +30);
   frame->SetXTitle(title_x); 
   frame->SetYTitle(title_y);
   frame->SetAxisRange(1E0,1E5,"Y");
  
   

   TH1D* hsum[nCan];
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
         
   TH1D *histo[nPlots];
   
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  
   for(Int_t i=0; i<nPlots; i++){
  Int_t ihist=(j*nPlots)+i;
  hname=Form("histo_%i",ihist);
  histo[i] = histo3d->ProjectionX(hname, histo3d->GetYaxis()->FindBin(ang1[i]), histo3d->GetYaxis()->FindBin(ang2[i]),histo3d->GetZaxis()->FindBin(rpang1[j]), histo3d->GetZaxis()->FindBin(rpang2[j]));
  //histo[i]->Sumw2();
  if(i==0) hsum[j]=(TH1D*)histo[0]->Clone("hsum");
  else hsum[j]->Add(histo[i]);
  //histo[i]->Scale(1/histo[i]->Integral());
   histo[i]->SetLineColor(i+1);
   histo[i]->SetLineWidth(2);
   histo[i]->DrawCopy("esame");
   legspectra->AddEntry(histo[i], Form("#phi-#Psi: %.2lf#pi-%.2lf#pi",ang1[i]/pi,ang2[i]/pi), "lp");
     
  }//i-loop
     legspectra->DrawClone("same");

  //TString output=Form("obr/%s/%s_%i.png",trigger.Data(),figName.Data(),j);
  //c1->SaveAs(output);
  
  if(j==0) histosum=(TH1D*)hsum[0]->Clone("hsum");
  else histosum->Add(hsum[j]);
   
  }//j-loop
  
  TCanvas *c2= new TCanvas("c2","",10,10,1500,900);
  c2->SetFillColor(0);
  c2->SetFrameFillColor(0);
  c2->SetGrid();
  c2->cd();
   gPad->SetLogy();
   TString htitle= "p_{T} vs #Psi_{EP}";
   frame->SetTitle(htitle);
   frame->DrawCopy("");
   
   TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  
   for(Int_t j=0; j<nCan; j++){
     hsum[j]->SetLineColor(j+1);
     hsum[j]->SetLineWidth(2);
   hsum[j]->DrawCopy("esame");
   legspectra->AddEntry(hsum[j], Form("#Psi: %.2lf#pi-%.2lf#pi",rpang1[j]/pi,rpang2[j]/pi), "lp");
   cout<<j<<" int "<<hsum[j]->Integral()<<endl;
   }
   legspectra->DrawClone("same");
   
 

   TCanvas *c3= new TCanvas("c3","",10,10,1500,900);
  c3->SetFillColor(0);
  c3->SetFrameFillColor(0);
  c3->SetGrid();
  c3->cd();
   gPad->SetLogy();
   htitle= Form("p_{T}^{corr}, p_{T}^{leading}>%i",pTlead);
   frame->SetTitle(htitle);
   frame->DrawCopy("");
   
   histoold1->SetLineColor(1);
   histoold2->SetLineColor(2);
   histoold1->DrawCopy("esame");
   
   //histoold2->DrawCopy("esame");
   histosum->SetLineColor(2);
   histosum->DrawCopy("esame");
   
   cout<<"integral old: "<<histoold2->Integral()<<" new: "<<histosum->Integral()<<endl;
   
      TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(histoold1,"old","lp");
  //legspectra->AddEntry(histoold2,"old2","lp");
  legspectra->AddEntry(histosum,"v2 reweighted","lp");
   legspectra->DrawClone("same");
   
  TString output=Form("obr/%s/%s_pTlead%i.png",trigger.Data(),figName.Data(),pTlead);
  c3->SaveAs(output);
  
   fout->cd();
   histosum->Write(Form("hmeasured_%i",pTlead));
   f->Close();
   fout->Close();
  }
     
   
