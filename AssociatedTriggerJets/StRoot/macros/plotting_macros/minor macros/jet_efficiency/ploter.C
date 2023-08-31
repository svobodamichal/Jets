// Example for fitting signal/background. 
// This example can be executed with:
// root > .x FittingDemo.C  (using the CINT interpreter)
// root > .x FittingDemo.C+ (using the native complier via ACLIC)
// function is : f(x) = (ax+b)exp(c(x-1.73)) + d*exp(-[(x-xo)/2sigma ]^2)   
//#include "TH1.h"
//#include "TH2.h"
//#include "TMath.h"
//#include "TF1.h"
//#include "TLegend.h"
//#include "TCanvas.h"
//#include "TFile.h"
//#include "TLatex.h"
//#include "TTree.h"
//#include "Riostream.h"
//#include "TStyle.h"

void plot(int pTlead=5) {
  //void plot(TString quant="mass",TString cut="", Int_t doFit=0,TString input="results_minimctree_single.root") {
    

 /* Int_t binX;
  Float_t infX;
  Float_t supX;*/
  TString title_x;
  Int_t binY;
  Float_t infY;
  Float_t supY;
  TString title_y;
  TString t1;
  TString t2;
  TString label;
  TString cut2;
  Double_t pi;
  Double_t norm;
  
  TString type="pythia - pythia corrected";

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptDate(1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleOffset(0.95,"Y");
  gStyle->SetTitleSize(0.07,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.03,"X");
  gStyle->SetLabelSize(0.03,"Y");
  
  
  //canvas size
  Int_t can_x=1200;
  Int_t can_y=800;
  
  Int_t maxX=100;
  
  pi=3.1415;

  Color_t colorList[30]={kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kOrange+2,kYellow+2,kBlue-4,kGreen,kOrange,kRed+2,kBlack-2,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.035);

  //TString input=Form("epsilon_R0.2_pTlead%i_20GeV_narrowmatrix.root",pTlead);
  //TString input=Form("epsilon_R0.2_pTlead%i_highpTjet.root",pTlead);
  //TString input=Form("epsilon_R0.2_pTlead%i_onejet.root",pTlead);
  //TString input=Form("epsilon_R0.2_pTlead%i_randomjet.root",pTlead);
  //TString input=Form("epsilon_R0.2_pTlead%i_ufrag.root",pTlead);
  //TString input=Form("epsilon_R0.2_pTlead%i_receffi_pTsmear.root",pTlead);
  TString input=Form("epsilon_R0.2_pTlead%i_effi_pTsmear.root",pTlead);
  TFile *f =new TFile(input);
  f->cd(); 
  
  
  TH1D* hpart=f->Get("hpart");
  TH1D* hpart2=f->Get("hpart2");
  TH1D* hdete=f->Get("hdete_measured");
  TH1D* hdp=f->Get("hdete_smeared");
  TH1D* hdu=f->Get("hdete_unfolded");
  TH1D* hepsilon=f->Get("hepsilon_transposed");
  TH1D* hepsilon2=f->Get("hepsilon_unfolded");
  TH2D* hresponse=f->Get("hresponse");
  TH2D* hunfold=f->Get("hunfoldM");
  TH2D* hunit=f->Get("hunit2");
  TH2D* htunit=f->Get("hunit");
  
  TH1D* projectX=(TH1D*)hresponse->ProjectionX();
  TH1D* projectY=(TH1D*)hresponse->ProjectionY();
  
  //TH2D* histo2[nprobes];
  //for(Int_t iprobe=0; iprobe<nprobes;iprobe++){
  
    TH1I* frame=new TH1I("frame","",1000,-100,100);

   TCanvas *c1 = new TCanvas("c1","",10,10,can_x,can_y);
   c1->cd();
   c1->SetFillColor(0);
   c1->SetFrameFillColor(0);
   c1->SetGrid();
   gPad->SetLogy();
   frame->SetAxisRange(0,maxX,"x");
   frame->SetMaximum(1E6);
   frame->SetMinimum(2E3);
   frame->GetXaxis()->SetTitle("jet p_{T}^{charged} (GeV/c)");
   frame->GetYaxis()->SetTitle("counts");
   frame->SetTitle(Form("R=0.2, p_{T}^{leading}>%i GeV/c",pTlead));

   hpart->SetLineColor(kRed);
   hdete->SetLineColor(kBlack);
   hdp->SetLineColor(kBlue);
   hdu->SetLineColor(kMagenta);
 
   hpart->SetLineWidth(2); 
   hdete->SetLineWidth(2); 
   hdp->SetLineWidth(2); 
   hdu->SetLineWidth(2); 

   frame->DrawCopy("e");
   hpart->DrawCopy("esame");
   //hpart2->DrawCopy("esame");
   hdete->DrawCopy("esame");
   hdu->DrawCopy("esame");
   hdp->DrawCopy("esame");

  /* 
   //projectX->Scale(hdete->Integral()/projectX->Integral());
   projectX->SetLineColor(kBlue);
   projectX->SetLineWidth(2);
   //projectY->Scale(hpart->Integral()/projectY->Integral());
   projectY->SetLineColor(kMagenta);
   projectY->SetLineWidth(2);
  
   projectX->Draw("same");
   projectY->Draw("same");*/

 /*
    TLegend *model_info = new TLegend(0.16037, 0.5425, 0.3848, 0.8981);
  model_info->SetFillStyle(0);
  model_info->SetBorderSize(0);
  model_info->SetMargin(0.05);
  model_info->SetHeader("Run11 AuAu 200 GeV/c, MB");
  model_info->AddEntry("", "charged jets", "");
  model_info->AddEntry("", "0-10% Central Collisions", "");
  model_info->AddEntry("", Form("N_{events} = %.1lfM", nEvents/1E6), "");
  model_info->AddEntry("", "Anti-k_{T} \t \t R = 0.3", "");
//model_info->AddEntry("", Form("p_{T}^{const} > 0.2 GeV/c", pTcut), "");
  model_info->AddEntry("", "A_{reco jet} > 0.2sr", "");
  //model_info->AddEntry("", "single particle fragmentation", "");
//model_info->AddEntry("", Form("p_{T}^{leading} > %.1lf GeV/c", pTthresh), "");
  model_info->DrawClone("same");
  */
  
  TLegend *leg = new TLegend(0.45, 0.64, 0.85, 0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.05);
  leg->SetTextSize(0.04);
  leg->AddEntry(hpart, Form("particle level, I=%.0lf.10^{3}",hpart->Integral("width")/1000.), "lp");
  //leg->AddEntry(hpart2, Form("particle level, I=%.0lf.10^{3}",hpart2->Integral("width")/1000.), "lp");
  leg->AddEntry(hdete, Form("detector level, I=%.0lf.10^{3}",hdete->Integral("width")/1000.), "lp");
  leg->AddEntry(hdu, Form("detector level - unfolded, I=%.0lf.10^{3}",hdu->Integral("width")/1000.), "lp");
  leg->AddEntry(hdp, Form("R^{T} x detector level, I=%.0lf.10^{3}",hdp->Integral("width")/1000.), "lp");
  /*leg->AddEntry(projectY, "RM - projection to Y (=particle level)", "lp");
  leg->AddEntry(projectX, "RM - projection to X (=detecor level)", "lp");*/

  leg->DrawClone("same");

  TString output=Form("obr/jeteffi_spectra_pTl%i.gif",pTlead);
  c1->SaveAs(output);

  TCanvas *c2 = new TCanvas("c2","",10,10,can_x,can_y);
  c2->cd();
  c2->SetFillColor(0);
  c2->SetFrameFillColor(0);
  c2->SetGrid();
  //gPad->SetLogy();
  frame->SetAxisRange(0,maxX,"x");
  frame->SetMaximum(2);
  frame->SetMinimum(0.2);
  frame->GetXaxis()->SetTitle("p_{T}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("efficiency");
  frame->SetTitle(Form("R=0.2, p_{T}^{leading}>%i GeV/c",pTlead));
  hepsilon->SetLineWidth(2);
  hepsilon2->SetLineWidth(2);
  frame->DrawCopy("");
  hepsilon2->SetLineColor(kRed);
  hepsilon->Draw("same");
  hepsilon2->Draw("same");

  TLegend *leg = new TLegend(0.45, 0.64, 0.85, 0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.05);
  leg->SetTextSize(0.05);
  leg->AddEntry(hepsilon, "R^{T} x p_{T}^{dete}/p_{T}^{part}", "lp");
  leg->AddEntry(hepsilon2, "unfolded p_{T}^{dete}/p_{T}^{part}", "lp");
  leg->DrawClone("same");
    
  TLine *one = new TLine(0, 1, maxX, 1);
   one->SetLineWidth(2);
    one->SetLineStyle(2);
   one->SetLineColor(kBlack);
  one->DrawClone("same");	

  output=Form("obr/epsilon_pTlead%i.gif",pTlead);
  c2->SaveAs(output);

  TCanvas *c3 = new TCanvas("c3","",10,10,can_x,can_y);
  c3->cd();
  c3->SetFillColor(0);
  c3->SetFrameFillColor(0);
  //c3->SetGrid();
  gPad->SetLogz();
  hresponse->SetTitle("Response Matrix");
  hresponse->GetXaxis()->SetTitle("p_{T} detector level");
  hresponse->GetYaxis()->SetTitle("p_{T} particle level");
  hresponse->SetAxisRange(0,maxX,"x");
  hresponse->SetAxisRange(0,maxX,"y");
  hresponse->Draw("COLZ");

  output=Form("obr/responsem_%i.gif",pTlead);
  c3->SaveAs(output);
  
  TCanvas *c4 = new TCanvas("c4","",10,10,can_x,can_y);
  c4->cd();
  c4->SetFillColor(0);
  c4->SetFrameFillColor(0);
  //c3->SetGrid();
  gPad->SetLogz();
  hunfold->SetTitle("Inverted Response Matrix");
  hunfold->GetXaxis()->SetTitle("p_{T} particle level");
  hunfold->GetYaxis()->SetTitle("p_{T} detector level");
  hunfold->SetAxisRange(0,maxX,"x");
  hunfold->SetAxisRange(0,maxX,"y");
  hunfold->Draw("COLZ");

  output=Form("obr/unfold_%i.gif",pTlead);
  c4->SaveAs(output);
  
  TCanvas *c5 = new TCanvas("c5","",10,10,can_x,can_y);
  c5->cd();
  c5->SetFillColor(0);
  c5->SetFrameFillColor(0);
  //c3->SetGrid();
  //gPad->SetLogz();
  hunit->SetTitle("RxR^{-1}");
  hunit->GetXaxis()->SetTitle("");
  hunit->GetYaxis()->SetTitle("");
  hunit->SetAxisRange(0,maxX,"x");
  hunit->SetAxisRange(0,maxX,"y");
  hunit->Draw("COLZ");

  output=Form("obr/unit_%i.gif",pTlead);
  c5->SaveAs(output);
  
  TCanvas *c6 = new TCanvas("c6","",10,10,can_x,can_y);
  c6->cd();
  c6->SetFillColor(0);
  c6->SetFrameFillColor(0);
  //c3->SetGrid();
  gPad->SetLogz();
  htunit->SetTitle("RxR^{T}");
  htunit->GetXaxis()->SetTitle("");
  htunit->GetYaxis()->SetTitle("");
  htunit->SetAxisRange(0,maxX,"x");
  htunit->SetAxisRange(0,maxX,"y");
  htunit->Draw("COLZ");

  output=Form("obr/unitT_%i.gif",pTlead);
  c6->SaveAs(output);
    
  
  TH1D* hdivp=hpart->Clone("hdivp");
  TH1D* hdivd=hdete->Clone("hdivd");
  //projectY->Rebin(2);
  //projectX->Rebin(2);
  projectY->SetAxisRange(0,maxX,"x");
  projectX->SetAxisRange(0,maxX,"x");
  TH1D* hpx=hdivp->Clone("hpx");
  TH1D* hpy=hdivd->Clone("hpy");
  hpx->Reset("MICE");
  hpy->Reset("MICE");
  for(int bin=1;bin<=hpx->GetNbinsX();bin++)
  {
    double pT=hpx->GetBinCenter(bin);
    int bn=projectX->FindBin(pT);
    double val=projectX->GetBinContent(bn);
    double err=projectX->GetBinError(bn);
    hpx->SetBinContent(bin,val);
    hpx->SetBinError(bin,err);
    
    pT=hpy->GetBinCenter(bin);
    bn=projectY->FindBin(pT);
    val=projectY->GetBinContent(bn);
    err=projectY->GetBinError(bn);
    hpy->SetBinContent(bin,val);
    hpy->SetBinError(bin,err);
 }
    
  hpy->Divide(hdivp);
  hpx->Divide(hdivd);

   TCanvas *c7 = new TCanvas("c7","",10,10,can_x,can_y);
   c7->cd();
   c7->SetFillColor(0);
   c7->SetFrameFillColor(0);
   c7->SetGrid();
   //gPad->SetLogy();
   frame->SetAxisRange(0,maxX,"x");
   frame->SetMaximum(1.2);
   frame->SetMinimum(0.4);
   frame->GetXaxis()->SetTitle("jet p_{T}^{charged} (GeV/c)");
   frame->GetYaxis()->SetTitle("ratio");
   frame->SetTitle(Form("R=0.2, p_{T}^{leading}>%i GeV/c",pTlead));

   hpy->SetLineColor(kRed);
   hpx->SetLineColor(kBlack);
 
   hpy->SetLineWidth(2); 
   hpx->SetLineWidth(2); 
    
   frame->DrawCopy("e");
   hpy->DrawCopy("esame");
   hpx->DrawCopy("esame");


  TLegend *leg = new TLegend(0.45, 0.25, 0.85, 0.5);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.05);
  leg->SetTextSize(0.04);
  leg->AddEntry(hpy, "RM projection to Y/particle level", "lp");
  leg->AddEntry(hpx, "RM projection to X/detector level", "lp");

  leg->DrawClone("same");

  TString output="obr/ratios.gif";
  c7->SaveAs(output);
  
     TCanvas *c8 = new TCanvas("c8","",10,10,can_x,can_y);
   c8->cd();
   c8->SetFillColor(0);
   c8->SetFrameFillColor(0);
   c8->SetGrid();
   gPad->SetLogy();
   frame->SetAxisRange(0,maxX,"x");
   frame->SetMaximum(1E6);
   frame->SetMinimum(2E3);
   frame->GetXaxis()->SetTitle("jet p_{T}^{charged} (GeV/c)");
   frame->GetYaxis()->SetTitle("counts");
   frame->SetTitle(Form("R=0.2, p_{T}^{leading}>%i GeV/c",pTlead));

   hpart->SetLineColor(kRed);
   hdete->SetLineColor(kBlack);
   hpart->SetLineWidth(2); 
   hdete->SetLineWidth(2); 

   frame->DrawCopy("e");
   hpart->DrawCopy("esame");
   hdete->DrawCopy("esame");
   //projectX->Scale(hdete->Integral()/projectX->Integral());
   projectX->SetLineColor(kBlue);
   projectX->SetLineWidth(2);
   //projectY->Scale(hpart->Integral()/projectY->Integral());
   projectY->SetLineColor(kMagenta);
   projectY->SetLineWidth(2);
  
   projectX->Draw("same");
   projectY->Draw("same");

  
  TLegend *leg = new TLegend(0.45, 0.64, 0.85, 0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.05);
  leg->SetTextSize(0.04);
  leg->AddEntry(hpart, "particle level", "lp");
  leg->AddEntry(hdete, "detector level", "lp");
  leg->AddEntry(projectY, "RM - projection to Y (=particle level)", "lp");
  leg->AddEntry(projectX, "RM - projection to X (=detecor level)", "lp");

  leg->DrawClone("same");

  TString output=Form("obr/RM_project_pTl%i.gif",pTlead);
  c8->SaveAs(output);
}
