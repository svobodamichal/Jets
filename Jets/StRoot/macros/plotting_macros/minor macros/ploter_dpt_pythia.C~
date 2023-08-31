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

void plot(int doToymodel=0, TString trigger="MB") {
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
  gStyle->SetTitleSize(5,"Y");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleSize(0.9,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.09,"X");
  gStyle->SetLabelSize(0.03,"Y");
  
  //canvas size
  Int_t can_x=1200;
  Int_t can_y=800;
  
  pi=3.1415;

  Color_t colorList[30]={kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kOrange+2,kYellow+2,kBlue-4,kGreen,kOrange,kRed+2,kBlack-2,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.035);

  //const Int_t nprobes=11;
  const Int_t nprobes=1;
  //const Float_t embPt[nprobes]={0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 30.0};
  const Float_t embPt[nprobes]={30.0};
  TString input=Form("root/%s/histos_embeddedjet_R0.4_pythia.root",trigger.Data());
  if(doToymodel) input="root/toymodel/histos_dpTarea_R0.4_pTcut0.2_1Mcharged.root";
  TFile *f =new TFile(input);
  f->cd(); 
  Int_t nEvents=1E6;
  if(!doToymodel)
  {
  TH1I* hevents= (TH1I*) f->Get("hevts");
  nEvents=hevents->GetEntries();
  }
  TH1D* histo[nprobes];
  TH1D* histo0[nprobes]; //histograms with 1 +/- 1 bins removed
  TH1D* histop[nprobes]; //pythia dpT
  TH1D* histopc[nprobes]; //tracking efficiency corrected pythia dpT
  //TH2D* histo2[nprobes];
  for(Int_t iprobe=0; iprobe<nprobes;iprobe++){
  TString inHistName=Form("hdpT%.1lf_cut3", embPt[iprobe]);
  TString inHistNamep=Form("hdpT%.1lf_pyth_cut3", embPt[iprobe]);
  TString inHistNamepc=Form("hdpT%.1lf_pythCorr_cut3", embPt[iprobe]);
  if(doToymodel)inHistName=Form("hdpT_Acut_pTemb%.1lf", embPt[iprobe]);
  //TString inHistName2D=Form("hdpT_area_pTemb%.1lf", embPt[iprobe]);
  TString outHistName=Form("histo_%i",iprobe);
  TString outHistNamep=Form("histop_%i",iprobe);
  TString outHistNamepc=Form("histopc_%i",iprobe);
  //TString outHistName2D=Form("histo2D_%i",iprobe);
  //TH2D* histo2D = (TH2D*) f->Get(inHistName2D);
  TH1D* histo1D = (TH1D*) f->Get(inHistName);
  histo1D->Sumw2();
  histo[iprobe] = (TH1D*) histo1D->Clone(outHistName);
  histo0[iprobe] = (TH1D*) histo1D->Clone(outHistName);
  
  TH1D* histo1Dp = (TH1D*) f->Get(inHistNamep);
  histo1D->Sumw2();
  histop[iprobe] = (TH1D*) histo1Dp->Clone(outHistNamep);
  
  TH1D* histo1Dpc = (TH1D*) f->Get(inHistNamepc);
  histo1D->Sumw2();
  histopc[iprobe] = (TH1D*) histo1Dpc->Clone(outHistNamepc);
  // removing 1 +/- 1
  /*for(Int_t bin = 1; bin <= 2000; bin++){
   if(histo0[iprobe]->GetBinContent(bin) == histo0[iprobe]->GetBinError(bin))
     { 
       histo0[iprobe]->SetBinContent(bin, 0);
       histo0[iprobe]->SetBinError(bin, 0);
     }
    }
*/
  
  //histo2[iprobe] = (TH2D*) histo2D->Clone(outHistName2D);
  }

   TCanvas *c1 = new TCanvas("c1","",10,10,can_x,can_y);
   c1->cd();
   c1->SetFillColor(0);
   c1->SetFrameFillColor(0);
   c1->SetGrid();
   gPad->SetLogy();
   //histo[0]->SetTitle("#eta_{emb}=flat");
   //histo[0]->SetTitle("#deltap_{T} distributions");
   histo[0]->SetTitle(Form("#deltap_{T} distributions, p_{T}^{probe}=%.1lf", embPt[0]));
for(Int_t iPR=0; iPR<nprobes; iPR++){
   histo[iPR]->SetAxisRange(-40,50,"x");
   histo[iPR]->SetLineColor(colorList[iPR]);
   norm=histo[iPR]->Integral();
   histo[iPR]->Scale(1/norm);
   histo[iPR]->SetMaximum(0.5);
   histo[iPR]->SetMinimum(0.000001);
   histo[iPR]->SetLineWidth(2);
        if(iPR==0){
    histo[iPR]->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
    histo[iPR]->GetYaxis()->SetTitle("probability");
    histo[iPR]->DrawCopy("e");
   // histo[iPR]->Fit("gaus");
  }
  else {histo[iPR]->DrawCopy("esame");}
  
  
  histop[iPR]->SetLineColor(colorList[iPR+1]);
  norm=histop[iPR]->Integral();
  histop[iPR]->Scale(1/norm);
  histop[iPR]->SetLineWidth(2);
  histop[iPR]->DrawCopy("esame");
  
  
  histopc[iPR]->SetLineColor(colorList[iPR+2]);
  norm=histopc[iPR]->Integral();
  histopc[iPR]->Scale(1/norm);
  histopc[iPR]->SetLineWidth(2);
  histopc[iPR]->DrawCopy("esame");
   
}
 latex->DrawLatex(0.42, 0.3,"STAR Preliminary");
 
    TLegend *model_info = new TLegend(0.16037, 0.5425, 0.3848, 0.8981);
  model_info->SetFillStyle(0);
  model_info->SetBorderSize(0);
  model_info->SetMargin(0.05);
  model_info->SetHeader(Form("Run11 AuAu 200 GeV/c, %s",trigger.Data()));
  if(doToymodel)model_info->SetHeader("Toymodel: RHIC kinematics");
  if(doToymodel)model_info->AddEntry("","single particle fragmentation","");
  model_info->AddEntry("", "charged jets", "");
  model_info->AddEntry("", "0-10% Central Collisions", "");
  model_info->AddEntry("", Form("N_{events} = %.1lfM", nEvents/1E6), "");
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
for(Int_t iPR=0; iPR<nprobes; iPR++){	
   //leg->AddEntry(histo[iPR], Form(%.1lf GeV/c", embPt[iPR]), "l");
   
   leg->AddEntry(histo[iPR], Form("Single particle, %.1lf GeV/c", embPt[iPR]), "l");
   leg->AddEntry(histop[iPR], Form("Pythia, %.1lf GeV/c", embPt[iPR]), "l");
   leg->AddEntry(histopc[iPR], Form("Pythia corr., %.1lf GeV/c", embPt[iPR]), "l");
}
   //leg->AddEntry("","#sigma =~ 3.6 GeV/c","");
   leg->SetTextSize(0.04);
   leg->DrawClone("same");

  //TString output="obr/deltapt_embPt_all.png";
  //TString output=Form("obr/%s/deltapt_data_etaEmbFlat.png",trigger.Data());
  TString output=Form("obr/%s/deltapt_pythia_%.1lf.png",trigger.Data(),embPt[0]);
  if(doToymodel) output="obr/deltapt_toymodel_etaEmbFlat.png";
  c1->SaveAs(output);
  /*
  TCanvas *c2 = new TCanvas("c2","",10,10,can_x,can_y);
   c2->SetFillColor(0);
   c2->SetFrameFillColor(0);
   c2->SetGrid();
   gPad->SetLogy();
   histo0[0]->SetYTitle("probability");
   histo0[0]->SetTitle("without 1+/-1 bins");
for(Int_t iPR=0; iPR<nprobes; iPR++){
   histo0[iPR]->SetAxisRange(-40,50,"x");
   histo0[iPR]->SetLineColor(colorList[iPR]);
   norm=histo0[iPR]->Integral();
   histo0[iPR]->Scale(1/norm);
   histo0[iPR]->SetMaximum(0.5);
   histo0[iPR]->SetMinimum(0.000001);
        if(iPR==0){
    histo0[iPR]->DrawCopy("e");
  }
   else {histo0[iPR]->DrawCopy("esame");}
}

  leg->DrawClone("same");
  model_info->DrawClone("same");
  TString output=Form("obr/%s/deltapt_data_etaEmbFlat_no1pm1bins.png",trigger.Data());
  c2->SaveAs(output);
*/
}
