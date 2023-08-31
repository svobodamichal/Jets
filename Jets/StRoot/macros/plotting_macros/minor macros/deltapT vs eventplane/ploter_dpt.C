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

void plot(int pTlead=0, int probe=0) {
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
  /*gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.9,"Y");
  gStyle->SetTitleOffset(0.5,"Y");
  gStyle->SetTitleSize(0.9,"X");
  gStyle->SetTitleOffset(0.95,"X");
  gStyle->SetLabelSize(0.09,"X");
  gStyle->SetLabelSize(0.04,"Y");
  */
    gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.09);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.045,"Y");
  gStyle->SetTitleOffset(0.95,"Y");
  gStyle->SetTitleSize(0.07,"X");
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

  const Int_t nprobes=16;
  //const Int_t nprobes=11;
  //const Float_t embPt[nprobes]={0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 30.0};
  
 const Float_t embPt[nprobes]= {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0};
  TString input="histos_embeddedjet_R0.3.root";
  TFile *f =new TFile(input);
  f->cd(); 
  
  TH1I* hevents= (TH1I*) f->Get("hevts");
  Int_t nEvents=hevents->GetEntries();
  
  TH1D* histoin;
  TH1D* histoout;
  //TH2D* histo2[nprobes];
  //for(Int_t iprobe=0; iprobe<nprobes;iprobe++){
  
 
  TString inHistNameIn=Form("delta_pt_inplane_%i", pTlead);
  TString inHistNameOut=Form("delta_pt_outplane_%i", pTlead);

  TString outHistNameIn=Form("histoin_%i",probe);
  TString outHistNameOut=Form("histoout_%i",probe);
  
  TH2D* histo2Din = (TH2D*) f->Get(inHistNameIn);
  histo2Din->Sumw2();
  Int_t bin= histo2Din->GetXaxis()->FindBin(embPt[probe]);
  histoin = (TH1D*) histo2Din->ProjectionY(outHistNameIn,bin,bin);
  histoin->Scale(1.0/histoin->Integral());
  
  TH2D* histo2Dout = (TH2D*) f->Get(inHistNameOut);
  histo2Dout->Sumw2();
  bin= histo2Dout->GetXaxis()->FindBin(embPt[probe]);
  histoout = (TH1D*) histo2Dout->ProjectionY(outHistNameOut,bin,bin);
  histoout->Scale(1.0/histoout->Integral());
  
  //}//


   TCanvas *c1 = new TCanvas("c1","",10,10,can_x,can_y);
   c1->cd();
   c1->SetFillColor(0);
   c1->SetFrameFillColor(0);
   c1->SetGrid();
   gPad->SetLogy();
   //histo[0]->SetTitle("#eta_{emb}=flat");
   histoout->SetAxisRange(-30,30,"x");
   histoout->SetMaximum(0.5);
   histoout->SetMinimum(1E-7);
   //histoout->SetTitleSize(0.2,"X");
   //histoout->SetLabelSize(0.09,"X");
   histoout->GetXaxis()->SetTitle("#deltap_{T} (GeV/c)");
   histoout->GetYaxis()->SetTitle("probability");
   histoout->SetTitle(Form("#deltap_{T} distributions, SP probe p_{T} = %.1lf GeV/c",embPt[probe]));

   histoin->SetLineColor(kRed);
   histoout->SetLineColor(kBlack);
 
   histoin->SetLineWidth(2); 
   histoout->SetLineWidth(2); 
  
   histoout->DrawCopy("e");
   histoin->DrawCopy("esame");
   
   TF1* fgaus=new TF1("fgaus","[0]*exp(-((x-[1])*(x-[1])/(2.*[2]*[2])))",-10,0);
   fgaus->SetParameters(0.01,0.0,4.0);

   fgaus->SetLineColor(kGreen);
   histoin->Fit("fgaus","R0");
   Double_t par[3];
   fgaus->GetParameters(par);
   Double_t muin=par[1];
   fgaus->DrawCopy("same");
   
   fgaus->SetLineColor(kBlue);
   histoout->Fit("fgaus","R0");
   fgaus->GetParameters(par);
   Double_t muout=par[1];
   fgaus->DrawCopy("same");
     

      
   Double_t deltam=muin-muout;
   cout<<"mean in: "<<muin<<"mean out: "<<muout<<" delta: "<<deltam<<endl;
 
 latex->DrawLatex(0.42, 0.3,"STAR Preliminary");
 
 
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
  
  
  TLegend *leg = new TLegend(0.70, 0.75, 0.95, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.05);
  leg->SetTextSize(0.04);
  leg->AddEntry(histoin, "in plane", "l");
  leg->AddEntry(histoout, "out of plane", "l");

  leg->DrawClone("same");

  TString output=Form("obr/deltapt_eventplane_%i.gif",probe);
  c1->SaveAs(output);
  
  TCanvas *c2 = new TCanvas("c2","",10,10,can_x,1.3*can_y);
  TPad*    upperPad = new TPad("upperPad", "upperPad", .005, 0.301, .995, .995);
  TPad*    lowerPad = new TPad("lowerPad", "lowerPad", 
			       .005, .005, .995, .3);
  upperPad->Draw(); 			       
      lowerPad->Draw();

   //c2->Divide(1,2);
   upperPad->cd();
   c2->SetFillColor(0);
   c2->SetFrameFillColor(0);
   c2->SetGrid();
   upperPad->SetLogy();
   upperPad->SetGrid();
   //gPad->SetCanvasSize(can_x,1.*can_y);
   TH1D* hshift=(TH1D*) histoin->Clone("hshift");
   hshift->Reset("MICE");
   Int_t nbins=hshift->GetNbinsX();
   for(Int_t bin=1; bin<=nbins; bin++){
   Double_t pTold=histoin->GetBinCenter(bin);
   Double_t yield=histoin->GetBinContent(bin);
   if(!yield>0)continue;
   Double_t error=histoin->GetBinError(bin);
   Double_t pTnew=pTold-deltam;
   Int_t binnew=hshift->FindBin(pTnew);
   hshift->SetBinContent(binnew,yield);
   hshift->SetBinError(binnew,error);
  }
  histoout->DrawCopy();
  hshift->DrawCopy("esame");
  histoout->DrawCopy("esame");
   
  model_info->DrawClone("same");
  
  TLegend *leg = new TLegend(0.70, 0.75, 0.95, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.05);
  leg->SetTextSize(0.04);
  leg->AddEntry(hshift, "in plane - shifted", "l");
  leg->AddEntry(histoout, "out of plane", "l");

  leg->DrawClone("same");
  
  //c2->cd(2)->SetCanvasSize(can_x,0.3*can_y);
  lowerPad->cd();

  //gPad->SetLogy();
  lowerPad->SetGrid();
  //gPad->SetCanvasSize(can_x,0.3*can_y);
  TH1D* hratio=(TH1D*) hshift->Clone();
  hratio->Divide(histoout);
  hratio->SetTitle("");
  hratio->GetYaxis()->SetTitleSize(0.085);
  hratio->GetYaxis()->SetLabelSize(0.05);
  hratio->GetYaxis()->SetTitleOffset(0.6);
  hratio->SetAxisRange(-30,30,"x");
  hratio->SetMaximum(3);
  hratio->SetMinimum(0.2);
  hratio->GetYaxis()->SetTitle("in/out");
  hratio->DrawCopy("e");
  
  TLine *one = new TLine(-30, 1, 30, 1);
  one->SetLineWidth(2);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  
  output=Form("obr/deltapt_eventplane_%i_shift.gif",probe);
  c2->SaveAs(output);

}
