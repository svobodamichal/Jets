#include <stdio.h>
#include "iomanip.h"
#include "TString.h"

void plotptleadratio2(TString version="GPC2",TString ext="pdf") {
  gROOT->Reset();
  gROOT->Clear();
  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleYSize(0.04);
  
	TString outputDir=Form("../../plotting_out/obr/%s/pp",version.Data());
    TString outputDir2="obr";
	
  TFile *f1 = new TFile("data/jethistos_370.root");

  TString cname,hname,hname7,hname5,hname0;
  TString cname_txt, hname_txt;
  int ptlead[4]={0,3,5,7};
  double rpar[5]={0.2,0.3,0.4,0.5,0.6};

  TCanvas *c=new TCanvas("c","pt leading spectra ratios",1200,400);
  c->Divide(3,1);

  
  for(int ir=0;ir<3;ir++) {
    c->cd(ir+1);
    gPad->SetLogy(1);
    
    hname="hchjet_pTleadrat_R0";
    hname+=ir+2;
    
    hname_txt="R=0.";
    hname_txt+=ir+2;
    
    TH1D *h1 = new TH1D(hname,hname_txt, 100,0,100); //with the fudge factor applied
    h1->Sumw2();
    
    hname7="hchjet_pT_R0";
    hname7+=ir+2;
    hname7+="_pTl7_sum_fudge";
    
    hname5="hchjet_pT_R0";
    hname5+=ir+2;
    hname5+="_pTl5_sum_fudge";
    
    
    h1->Divide(((TH1D*)f1->Get(hname7)),((TH1D*)f1->Get(hname5)),1.,1.,"e");
    h1->Rebin(2);
    h1->Scale(0.5);
   
    h1->SetXTitle("p_{T,jet}^{charged} (GeV/c)");
    h1->SetYTitle("p_{T,lead}^{min} > 7 GeV/c / p_{T,lead}^{min} > 5 GeV/c");
    h1->GetXaxis()->SetRangeUser(0.,50.);
    h1->SetMarkerStyle(20);
    h1->Draw("CL");
  }
  c->Print(Form("%s/ptleading_ratios_pythia.%s",outputDir2.Data(),ext.Data()));



  
  TCanvas *c1=new TCanvas("c1","pt leading spectra ratios",500,400);
  c1->cd();

  TLegend *leg = new TLegend(0.35,0.2,0.5,0.45,"p+p #\sqrt{s} = 200 GeV charged jets spectra ratio");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  
  for(int ir=0;ir<3;ir++) {
    c1->cd(ir+1);
    gPad->SetLogy(1);
    
    hname="hchjet_pTleadrat_R0";
    hname+=ir+2;
    
    hname_txt="R=0.";
    hname_txt+=ir+2;
    
    TH1D *h1 = new TH1D(hname,hname_txt, 100,0,100); //with the fudge factor applied
    h1->Sumw2();
    
    hname7="hchjet_pT_R0";
    hname7+=ir+2;
    hname7+="_pTl7_sum_fudge";
    
    hname5="hchjet_pT_R0";
    hname5+=ir+2;
    hname5+="_pTl5_sum_fudge";
    
    
    h1->Divide(((TH1D*)f1->Get(hname7)),((TH1D*)f1->Get(hname5)),1.,1.,"e");
    h1->Rebin(2);
    h1->Scale(0.5);
   
    h1->SetXTitle("p_{T,jet}^{charged} (GeV/c)");
    h1->SetYTitle("p_{T,lead}^{min} > 7 GeV/c / p_{T,lead}^{min} > 5 GeV/c");
    h1->GetXaxis()->SetRangeUser(0.,50.);
    h1->GetYaxis()->SetRangeUser(0.1,1.5);
    
    h1->SetMarkerStyle(20+ir);
  
    if(ir==0) {
      h1->SetMarkerColor(1);
    }
    
    if(ir==1) {
      h1->SetMarkerColor(2);
    }
    
    if(ir==2) {
      h1->SetMarkerColor(4);
    }
    
    h1->Draw("same");
    
    leg->AddEntry(h1,hname_txt,"P");
    
  }
 
  
  leg->Draw();
 
  c1->Print(Form("%s/ptleading_ratios_pythia_one.%s",outputDir2.Data(),ext.Data()));


  //plot ratio 5/0
  TCanvas *c2=new TCanvas("c2","pt leading spectra ratios",500,400);
  c2->cd();

  
 
  TLegend *leg = new TLegend(0.42,0.15,0.8,0.4,"    p+p #\sqrt{s} = 200 GeV");
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetNColumns(3);
 
  // leg->AddEntry((TObject*)0," "," ");
  leg->AddEntry((TObject*)0,"                 p_{T,cut1}/p_{T_cut2} (GeV/#it{c})"," ");
  leg->AddEntry((TObject*)0," "," ");
  leg->AddEntry((TObject*)0," "," ");
  
  gPad->SetLogy(1);
  TH2F *frame = new TH2F("ratio","",35, 0.,35.,10,0.3,1.5);
  frame->GetYaxis()->SetTitleSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.035);
  frame->GetXaxis()->SetTitleOffset(1.2); //Jana
  frame->SetXTitle("p_{T,jet}^{charged} (GeV/#it{c})");
  frame->SetYTitle("dN_{ch,jet}/dp_{T}d#\eta (p_{T,lead}^{min}>p_{T,cut1})/dN_{ch,jet}/dp_{T}d#\eta (p_{T,lead}^{min}>p_{T,cut2})");
  frame->Draw();

    
  for(int ir=0;ir<3;ir++) {
  
    
    hname="hchjet_pTleadrat_R0";
    hname+=ir+2;
    
    hname_txt="R=0.";
    hname_txt+=ir+2;
    
    TH1D *h1 = new TH1D(hname,"", 100,0,100); //with the fudge factor applied 7/5
    h1->Sumw2();

    TH1D *h2 = new TH1D(hname,"", 100,0,100); //with the fudge factor applied 5/0
    h2->Sumw2();

    hname7="hchjet_pT_R0";
    hname7+=ir+2;
    hname7+="_pTl7_sum_fudge";
    
    hname5="hchjet_pT_R0";
    hname5+=ir+2;
    hname5+="_pTl5_sum_fudge";
    
    hname0="hchjet_pT_R0";
    hname0+=ir+2;
    hname0+="_pTl0_sum_fudge";
    
    
    h1->Divide(((TH1D*)f1->Get(hname7)),((TH1D*)f1->Get(hname5)),1.,1.,"e");
    h1->Rebin(2);
    h1->Scale(0.5);
    
    
    // h1->SetXTitle("p_{T,jet}^{charged} (GeV/c)");
    //  h1->SetYTitle("p_{T,lead}^{min} > 5 GeV/c / p_{T,lead}^{min} > 0 GeV/c");
    // h1->SetYTitle("p_{T,lead}^{min} > p_{T,cut1} GeV/c / p_{T,lead}^{min} > p_{T,cut2} GeV/c");
    h1->GetXaxis()->SetRangeUser(0.,35.);
    h1->GetYaxis()->SetRangeUser(0.3,1.5);
    


    h2->Divide(((TH1D*)f1->Get(hname5)),((TH1D*)f1->Get(hname0)),1.,1.,"e");
    h2->Rebin(2);
    h2->Scale(0.5);
   
    // h2->SetXTitle("p_{T,jet}^{charged} (GeV/c)");
    // h2->SetYTitle("p_{T,lead}^{min} > p_{T,cut1} GeV/c / p_{T,lead}^{min} > p_{T,cut2} GeV/c");
    h2->GetXaxis()->SetRangeUser(0.,35.);
    h2->GetYaxis()->SetRangeUser(0.3,1.5);
    
  
    if(ir==0) {
      // h1->SetMarkerColor(2); //was 1
      //  h2->SetMarkerColor(2); //was 1
      //  h1->SetMarkerStyle(20);
      //  h2->SetMarkerStyle(24);
      h1->SetLineColor(2);
      h2->SetLineColor(2);
      
    }
    
    if(ir==1) {
      // h1->SetMarkerColor(2);
      // h2->SetMarkerColor(2);
      // h1->SetMarkerStyle(21);
      // h2->SetMarkerStyle(25);
      
    }
    
    if(ir==2) {
      // // h1->SetMarkerColor(4);
      // // h2->SetMarkerColor(4);
      // // h1->SetMarkerStyle(22);
      // // h2->SetMarkerStyle(26);

      // h1->SetMarkerColor(4);
      // h2->SetMarkerColor(4);
      // h1->SetMarkerStyle(21);
      // h2->SetMarkerStyle(25);
      
    }


    

    // if((ir==0)||(ir==2)) {
    // if(ir==2) {
    //  h1->Draw("c hist same");
    // h2->Draw("c hist same");
    //  h1->Fit("pol3"," " , " ", 6.,37.);
    //      h1->Draw("same");
    // h2->Fit("pol3"," ","",6.,30.);
    

    //double p1[4]={-0.599324,0.161885,-0.00562325,6.6387e-05}; 
    if(ir==0) {
      
      TF1 *fconst1 = new TF1("fconst1","[0]*(1.0-TMath::Exp([1]*x+[2]))",6,35);
      fconst1->SetRange(6.,35.);
      fconst1->SetParameters(1.0,-0.3,1.0);
      h1->Fit("fconst1","0"," ", 6.,35.);
      TF1 *fun1 = (TF1*)h1->GetListOfFunctions()->FindObject("fconst1");
      fun1->SetRange(5,35);
      fun1->SetLineColor(2);
      fun1->SetLineStyle(1);
      fun1->Draw("same");
      

      
      TF1 *fconst2 = new TF1("fconst2","[0]*(1.0-TMath::Exp([1]*x+[2]))",6,35);
      fconst2->SetRange(6.,35.);
      fconst2->SetParameters(1.0,-0.3,1.0);
      h2->Fit("fconst2","0"," ", 6.,35.);
      TF1 *fun2 = (TF1*)h2->GetListOfFunctions()->FindObject("fconst2");
      fun2->SetRange(5,35);
      fun2->SetLineColor(2);
      fun2->SetLineStyle(2);
      fun2->Draw("same");
      
      
      hname_txt="R = 0.";
      hname_txt+=ir+2;
      hname_txt+="         ";
      leg->AddEntry((TObject*)0,hname_txt," ");
      
      leg->AddEntry(fun1," 7/5                 ","L");
      leg->AddEntry(fun2," 5/0","L");
      
    }
      
     if(ir==2) {
      
      TF1 *fconst3 = new TF1("fconst3","[0]*(1.0-TMath::Exp([1]*x+[2]))",6,35);
      fconst3->SetRange(6.,35.);
      fconst3->SetParameters(1.0,-0.3,1.0);
      h1->Fit("fconst3","0"," ", 6.,35.);
      TF1 *fun3 = (TF1*)h1->GetListOfFunctions()->FindObject("fconst3");
      fun3->SetRange(5,35);
      fun3->SetLineColor(4);
      fun3->SetLineStyle(1);
      fun3->Draw("same");
      

      
      TF1 *fconst4 = new TF1("fconst4","[0]*(1.0-TMath::Exp([1]*x+[2]))",8,35);
      fconst4->SetRange(8.,35.);
      fconst4->SetParameters(1.0,-0.1,1.0);
      h2->Fit("fconst4","0"," ", 8.,35.);
      TF1 *fun4 = (TF1*)h2->GetListOfFunctions()->FindObject("fconst4");
      fun4->SetRange(5,35);
      fun4->SetLineColor(4);
      fun4->SetLineStyle(2);
      fun4->Draw("same");
      
      
      hname_txt="R = 0.";
      hname_txt+=ir+2;
      hname_txt+="         ";
      leg->AddEntry((TObject*)0,hname_txt," ");
      
      leg->AddEntry(fun3," 7/5                 ","L");
      leg->AddEntry(fun4," 5/0","L");
      
    }

   
  }
 
  
  leg->Draw();
  TLatex *t1 = new TLatex(16.6,0.6,"PYTHIA 6.428");
  t1->SetTextSize(0.045);
  t1->Draw("same");
  
  c2->Print(Form("%s/ptleading_ratios_pythia_combined_lines.%s",outputDir.Data(),ext.Data()));
 
}

