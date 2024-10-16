//comparison of 2 spectra

void compare(Float_t pTthresh=0.0)
{
  
  Float_t R=0.4;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
    //canvas size
  Int_t can_x=1050;
  Int_t can_y=700;
  
  Float_t iscale=1000; //integral scaler - so it can fit the legend window
  TString iscale_str="1E-3";
  
  TString str;

  TString wrkdir_true = "./root/toymodel/jetonly";

  //wrkdir_true = "./root/toymodel/jetonly";
  TString outdir = "./obr/";
 
    
  //const Int_t nProbes=7; //number of embedded probe types
  //Color_t color[] = {kGreen+3, kMagenta, kCyan+1, kBlue, kYellow+2};

  TLatex *latex = new TLatex();
  latex->SetNDC();
  Color_t colorList[8]={kBlack,kRed,kGreen+3,kBlue,kMagenta+2,kOrange+2,kYellow+3,kBlue+4};

  TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
 
  Int_t nevents=1E6; //number of events in toymodel simulation
  Double_t scale_true=1./((2-2*R)*2.*TMath::Pi()*nevents);

  //T_AA scaled pp spectrum
  str = Form("%s/histos_jets_charged_R%.1lf_pTcut0.2.root", wrkdir_true.Data(),R);
  //if(doToymodel==0) str = Form("%s/histos_jets_R%.1lf_pTcut0.2_10Mcharged.root", wrkdir_true.Data(),R);
  TFile *ftruth1 = new TFile(str.Data(), "OPEN");
  TH2D *hPtRecpTleadingPrior = (TH2D*)ftruth1->Get("fhPtRecpTleading");
  Int_t firstbin = hPtRecpTleadingPrior->GetXaxis()->FindBin(pTthresh);
  Int_t lastbin = hPtRecpTleadingPrior->GetNbinsX();
  TH1D *htemp1 = hPtRecpTleadingPrior->ProjectionY("htemp1", firstbin, lastbin);
  TH1D *htruth1 = (TH1D*)htemp1->Clone("htruth1");
  htruth1->Sumw2();
  htruth1->Reset("MICE");
  for(Int_t bin = 1; bin <= htemp1->GetNbinsX(); bin++)
    {
      Double_t pTtmp = htemp1->GetBinCenter(bin);
      Double_t yield = htemp1->GetBinContent(bin);
      //Double_t error = htemp->GetBinError(bin);

      htruth1->Fill(pTtmp, yield);
    }
  for(Int_t bin = 1; bin <= htruth1->GetNbinsX(); bin++)
    {
      Double_t yield = htruth1->GetBinContent(bin);
      htruth1->SetBinError(bin,TMath::Sqrt(yield));
    }
  htruth1->SetMarkerColor(kBlack);
  htruth1->SetLineColor(kBlack);
  htruth1->SetMarkerStyle(kOpenCircle);
  htruth1->Scale(scale_true, "width");
  htruth1->SetMarkerSize(0.6);
  
 str = Form("%s/histos_jets_full_R%.1lf_pTcut0.2.root", wrkdir_true.Data(),R);
  //if(doToymodel==0) str = Form("%s/histos_jets_R%.1lf_pTcut0.2_10Mcharged.root", wrkdir_true.Data(),R);
  TFile *ftruth2 = new TFile(str.Data(), "OPEN");
  hPtRecpTleadingPrior = (TH2D*)ftruth2->Get("fhPtRecpTleading");
  firstbin = hPtRecpTleadingPrior->GetXaxis()->FindBin(pTthresh);
  lastbin = hPtRecpTleadingPrior->GetNbinsX();
  TH1D *htemp2 = hPtRecpTleadingPrior->ProjectionY("htemp2", firstbin, lastbin);
  TH1D *htruth2 = (TH1D*)htemp2->Clone("htruth2");
  htruth2->Sumw2();
  htruth2->Reset("MICE");
  for(Int_t bin = 1; bin <= htemp2->GetNbinsX(); bin++)
    {
      Double_t pTtmp = htemp2->GetBinCenter(bin);
      Double_t yield = htemp2->GetBinContent(bin);
      //Double_t error = htemp->GetBinError(bin);

      htruth2->Fill(pTtmp, yield);
    }
  for(Int_t bin = 1; bin <= htruth2->GetNbinsX(); bin++)
    {
      Double_t yield = htruth2->GetBinContent(bin);
      htruth2->SetBinError(bin,TMath::Sqrt(yield));
    }
  htruth2->SetMarkerColor(kRed);
  htruth2->SetLineColor(kRed);
  htruth2->SetMarkerStyle(kOpenCircle);
  htruth2->Scale(scale_true, "width");
  htruth2->SetMarkerSize(0.6);
  
  

  // DRAWING SPECTRA
  TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
  cspectra->cd();
  cspectra->SetGrid();
  cspectra->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
  frame->GetYaxis()->SetTitle("1/N_{events} 1/2#pi d^{2}N/dp_{T}d#eta (GeV/c)^{-1}");
  frame->GetXaxis()->SetRangeUser(-30, 40);
  frame->GetYaxis()->SetRangeUser(1E-7, 1E1);
  str=Form("T_{AA}*d#sigma_{pp}/dp_{T}, p_{T}^{leading} > %.1lf GeV/c",pTthresh);
  frame->SetTitle(str);
  frame->DrawCopy("");

  htruth1->DrawCopy("E same");
  htruth2->DrawCopy("E same");

  TLegend *model_info = new TLegend(0.12037, 0.5625, 0.3348, 0.900);
  model_info->SetFillStyle(0);
  model_info->SetBorderSize(0);
  model_info->SetMargin(0.05);
  model_info->SetHeader("Toymodel: RHIC kinematics");
  //model_info->AddEntry("", "Single Particle", "");
  model_info->AddEntry("", "0-5% Central Collisions", "");
  model_info->AddEntry("", Form("N_{events} = %.1lfM", nevents/1E6), "");
  model_info->AddEntry("", Form("Anti-k_{T} \t \t R = %.1lf",R), "");
  model_info->AddEntry("", "p_{T}^{const} > 0.2 GeV/c", "");
  model_info->AddEntry("", "A_{reco jet} > 0.4sr", "");
  //model_info->AddEntry("", Form("p_{T}^{leading} > %.1lf GeV/c", pTthresh), "");
  //model_info->AddEntry("", Form("prior distribution: %s", prior_type[priorNo].Data()), "");
  model_info->DrawClone("same");

  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(htruth1,"charged jets", "lp");
  legspectra->AddEntry(htruth2,"full jets", "lp");
  legspectra->DrawClone("same");
  str = Form("%s/compare_true_%.1lf.png", outdir.Data(), pTthresh);
  cspectra->SaveAs(str.Data());
  



  

}
