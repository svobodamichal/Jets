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

void plot(TString type="BG_sp", TString type_long="Single particle fragmentation", TString trigger="MB", Double_t R=0.4) {
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
  
  

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptDate(1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.055,"Y");
  gStyle->SetTitleOffset(0.95,"Y");
  gStyle->SetTitleSize(0.06,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.03,"X");
  gStyle->SetLabelSize(0.03,"Y");
  
  //canvas size
  Int_t can_x=1200;
  Int_t can_y=800;
  
  pi=3.1415;

  Color_t colorList[30]={kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kOrange+2,kYellow+2,kBlue-4,kGreen,kOrange,kRed+2,kBlack-2,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.035);

  
  const Int_t nemb=10;
  //const Int_t nemb=11;
  //const Float_t fEmbPt[nemb]={0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 30.0};
 //const Double_t fEmbPt[] = {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 23.0, 26.0 ,29.0 , 32.0, 36.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 80.0, 90.0};
 const Double_t fEmbPt[] = {0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0 , 40.0, 70.0};
  TString input=Form("root/%s/R%.1lf/histos_embeddedjet_R%.1lf.root",trigger.Data(),R,R);
  TFile* finput = new TFile(input, "OPEN");
//cout<<"input set"<<endl;

  TString name = Form("delta_pt_%s",type.Data());
//cout<<"delta_pt hsmeargram name: "<<name<<endl;
  TH2D* hdpT2D = (TH2D*)finput->Get(name.Data());
//cout<<"delta_pt hsmeargram: "<<hdpT2D->GetTitle()<<endl;

      /*// removing 1 +/- 1
    Int_t nbins=400;
      for(Int_t bin = 1; bin <= nbins*nbins; bin++)
   if(hdpT2D->GetBinContent(bin) == hdpT2D->GetBinError(bin))
     {
       hdpT2D->SetBinContent(bin, 0);
       hdpT2D->SetBinError(bin, 0);
     }
*/
TH1D* hsmear[nemb];

   for(Int_t i=0; i<nemb; i++){
      Int_t bin = hdpT2D->GetXaxis()->FindBin(fEmbPt[i]);
      TString name=Form("hsmear_%i",i);
      hsmear[i]=hdpT2D->ProjectionY(name, bin, bin);
   }

TH1 *frame = new TH1I("frame", "", 1000, -100, +100);

   TCanvas *c1 = new TCanvas("c1","",10,10,can_x,can_y);
   c1->cd();
   c1->SetFillColor(0);
   c1->SetFrameFillColor(0);
   c1->SetGrid();
   gPad->SetLogy();
   //hsmear[0]->SetTitle("#eta_{emb}=flat");
   frame->SetTitle(Form("#deltap_{T}^{charged} distributions, %s", type_long.Data()));
   frame->SetAxisRange(-20,30,"x");
   frame->GetYaxis()->SetRangeUser(1E-5,1);
   frame->GetXaxis()->SetTitle("#deltap_{T}^{charged} (GeV/c)");
   frame->GetYaxis()->SetTitle("probability");
      frame->DrawCopy("");
for(Int_t iPR=0; iPR<nemb; iPR++){
   
   hsmear[iPR]->SetLineColor(colorList[iPR]);
   norm=hsmear[iPR]->Integral();
   hsmear[iPR]->Scale(1/norm);
   hsmear[iPR]->SetMaximum(1);
   hsmear[iPR]->SetMinimum(0.000001);
   hsmear[iPR]->SetLineWidth(2);
   hsmear[iPR]->DrawCopy("esame");
  
  /*
  hsmearp[iPR]->SetLineColor(colorList[iPR+1]);
  norm=hsmearp[iPR]->Integral();
  hsmearp[iPR]->Scale(1/norm);
  hsmearp[iPR]->SetLineWidth(2);
  hsmearp[iPR]->DrawCopy("esame");
  
  
  hsmearpc[iPR]->SetLineColor(colorList[iPR+2]);
  norm=hsmearpc[iPR]->Integral();
  hsmearpc[iPR]->Scale(1/norm);
  hsmearpc[iPR]->SetLineWidth(2);
  hsmearpc[iPR]->DrawCopy("esame");*/
   
}
 latex->DrawLatex(0.42, 0.3,"STAR Preliminary");
 
    TLegend *model_info = new TLegend(0.16037, 0.5425, 0.3848, 0.8981);
  model_info->SetFillStyle(0);
  model_info->SetBorderSize(0);
  model_info->SetMargin(0.05);
  model_info->SetHeader(Form("Run11 AuAu 200 GeV/c, %s",trigger.Data()));
  model_info->AddEntry("", "charged jets", "");
  model_info->AddEntry("", "0-10% Central Collisions", "");
  //model_info->AddEntry("", Form("N_{events} = %.1lfM", nEvents/1E6), "");
  model_info->AddEntry("", "Anti-k_{T} \t \t R = 0.4", "");
//model_info->AddEntry("", Form("p_{T}^{const} > 0.2 GeV/c", pTcut), "");
  model_info->AddEntry("", "A_{reco jet} > 0.4sr", "");
  //model_info->AddEntry("", "single particle fragmentation", "");
//model_info->AddEntry("", Form("p_{T}^{leading} > %.1lf GeV/c", pTthresh), "");
  model_info->DrawClone("same");
 
  TLegend *leg = new TLegend(0.70, 0.5425, 0.95, 0.8981);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetMargin(0.05);
   leg->AddEntry("", "probe p_{T}:", "");
for(Int_t iPR=0; iPR<nemb; iPR++){	
   leg->AddEntry(hsmear[iPR], Form("%.1lf GeV/c", fEmbPt[iPR]), "l");
   /*
   leg->AddEntry(hsmear[iPR], Form("Single particle, %.1lf GeV/c", fEmbPt[iPR]), "l");
   leg->AddEntry(hsmearp[iPR], Form("Pythia, %.1lf GeV/c", fEmbPt[iPR]), "l");
   leg->AddEntry(hsmearpc[iPR], Form("Pythia corr., %.1lf GeV/c", fEmbPt[iPR]), "l");*/
}
   //leg->AddEntry("","#sigma =~ 3.6 GeV/c","");
   leg->SetTextSize(0.04);
   leg->DrawClone("same");

  //TString output="obr/deltapt_fEmbPt_all.png";
  TString output=Form("obr/%s/deltapT/R%.1lf/deltapt_%s.png",trigger.Data(),R,type.Data());
  //TString output=Form("obr/%s/deltapt_pythia_%.1lf.png",trigger.Data(),fEmbPt[0]);
  c1->SaveAs(output);
  
}
