plot_priors(Float_t pTthresh=5.0,TString suff="")
{ 
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetOptDate(0);
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
  
    Int_t markerU=21; //marker style for unfolding
  Int_t markerP=20; //marker style for prior
  Int_t markerM=29; //marker style for measured
  Int_t markerB=22; //marker style for backfolded
  Int_t markerU2=25; //marker style for unfolding2
  Int_t markerP2=24; //marker style for prior2
  Int_t markerM2=30; //marker style for measured2
  Int_t markerB2=26; //marker style for backfolded2
  
  Float_t marker_size=1.2;
  Float_t line_width=2;
  
  Color_t colorList[]={kBlack,kRed,kGreen+3,kMagenta+2,kBlue,kOrange+2,kYellow+2,kBlue-4,kOrange-2,kGreen,kRed+2,7,kGreen+4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,28};
  
  //canvas size
  Int_t can_x=1200; //1600
  Int_t can_y=680; //900
  
  Double_t R[3]={0.3,0.4,0.2};
  TString prior_type[]={"flat","pythiadete","pythia","powlaw3","powlaw4","powlaw5","powlaw6","levy","levy_alex"};
  TString prior_type_name[]={"flat", "Pythia dete.","biased Pythia","1/pT^{3}","1/pT^{4}","1/pT^{5}","1/(pT)^{6}", "levy I", "levy"};
  Int_t priors[]={5,3,8,1,2};
  const Int_t npriors=5;
  Int_t npriors_py=0;
  
    Float_t spectraXmin=0;
  Float_t spectraXmax=40;
  Float_t spectraYmin=1E-7;
  Float_t spectraYmax=1E2;
  
  TString str;
  TFile *fpyt[3];
  TH1D *hprior_py[3];
  TH1D *hprior_py0[3];
  TH1D *hprior[npriors];
  
  for(int i=0;i<npriors_py;i++){
	str= Form("./root/histos_pythiajet_R%.1lf.root",R[i]);
	fpyt[i]= new TFile(str.Data(), "OPEN");
	TH2D *hPrior2d = (TH2D*)fpyt[i]->Get("hpT_pTlead");
	Int_t firstbin = hPrior2d->GetYaxis()->FindBin(pTthresh);
	Int_t zerobin = hPrior2d->GetYaxis()->FindBin(0);
	Int_t lastbin = hPrior2d->GetNbinsY();
	hprior_py[i] = (TH1D*)hPrior2d->ProjectionX(Form("prior_pythia_%i",i), firstbin, lastbin);
	hprior_py0[i] = (TH1D*)hPrior2d->ProjectionX(Form("prior_pythia0_%i",i), zerobin, lastbin);
	delete hPrior2d;
	hprior_py[i]->SetLineColor(colorList[i]);
    hprior_py[i]->SetMarkerColor(colorList[i]);
    hprior_py[i]->SetMarkerStyle(markerU);
    hprior_py[i]->SetMarkerSize(marker_size);
    hprior_py[i]->SetLineWidth(line_width);
	 hprior_py[i]->Rebin(4);
	 hprior_py[i]->Scale(1.0/hprior_py[i]->Integral(),"width");
	 hprior_py0[i]->Rebin(hprior_py0[i]->GetNbinsX()/hprior_py[i]->GetNbinsX());
	 hprior_py0[i]->Scale(hprior_py[i]->Integral(hprior_py[i]->FindBin(35),hprior_py[i]->GetNbinsX(),"width")/hprior_py0[i]->Integral(hprior_py0[i]->FindBin(35),hprior_py0[i]->GetNbinsX(),"width"));
  }
  
  //str = "./root/histos_prior_flatstart.root";
  str = "./root/histos_prior.root";

  TFile *fprior = new TFile(str.Data(), "OPEN");
  for(int i=0;i<npriors;i++){
	  str=Form("hprior_%s",prior_type[priors[i]].Data());
	  if(priors[i]==1 || priors[i]==2) str=Form("hprior_%s_R%.1lf",prior_type[priors[i]].Data(),R[1]);
	  TH2D *hprior2d = (TH2D*)fprior->Get(str.Data());
	  Int_t firstbin = hprior2d->GetYaxis()->FindBin(pTthresh);
	  Int_t lastbin = firstbin;
	  TString priorName=Form("prior_%i",priors[i]);
	  hprior[i] = hprior2d->ProjectionX(priorName,firstbin,lastbin,"e");
	  delete hprior2d;
	  hprior[i]->SetLineColor(colorList[i+npriors_py]);
      hprior[i]->SetMarkerColor(colorList[i+npriors_py]);
      hprior[i]->SetMarkerStyle(markerU);
      hprior[i]->SetMarkerSize(marker_size);
      hprior[i]->SetLineWidth(line_width);
		hprior[i]->Rebin(4);
	   hprior[i]->Scale(1.0/hprior[i]->Integral(),"width");
	  
  }
    
    TH1 *frame = new TH1I("frame", "", 1000, -100, +100);

  TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
  cspectra->cd();
  cspectra->SetGrid();
  cspectra->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("probability");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  //frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
  str="Priors";
  frame->SetTitle(str);
  frame->DrawCopy("");
  for(int i=0;i<npriors_py;i++){
	hprior_py[i]->DrawCopy("E same");
  }
  for(int i=0;i<npriors;i++){
	hprior[i]->DrawCopy("E same");
  }
  
  
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  for(int i=0;i<npriors_py;i++){
  legspectra->AddEntry(hprior_py[i], Form("pythia w/ p_{T}^{leading} cut (R=%.1lf)",R[i]), "lp");
  }
  for(int i=0;i<npriors;i++){
  legspectra->AddEntry(hprior[i], Form("%s",prior_type_name[priors[i]].Data()), "lp");
  }
  legspectra->DrawClone("same");
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.4, 0.8,Form("p_{T}^{leading}>%.1lf GeV/c",pTthresh));
    
  str = Form("priors_pTlead%.0lf_%s.gif",pTthresh,suff.Data());
  cspectra->SaveAs(str.Data());
 
 
  
  TH1D* hdiv_py[3];
  TH1D* hdiv_py0[3];
  TH1D* hdiv[npriors];
  for(int i=0;i<npriors_py;i++){
  hdiv_py[i]=(TH1D*)hprior_py[i]->Clone(Form("hdiv_py_%i",i));
  hdiv_py[i]->Divide(hprior[0]);
  hdiv_py0[i]=(TH1D*)hprior_py[i]->Clone(Form("hdiv_py0_%i",i));
  hdiv_py0[i]->Divide(hprior_py0[i]);
  }
  for(int i=1;i<npriors;i++){
  hdiv[i]=(TH1D*)hprior[i]->Clone(Form("hdiv_%i",i));
  hdiv[i]->Divide(hprior[0]);
  }
  
  TCanvas *cratio = new TCanvas("cratio","cratio",10,10,can_x,can_y);
  cratio->cd();
  cratio->SetGrid();
  cratio->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("ratio");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(1E-3, 1E4);
  str=Form("Priors - ratios wrt %s",prior_type_name[priors[0]].Data());
  frame->SetTitle(str);
  frame->DrawCopy("");
  for(int i=0;i<npriors_py;i++){
	hdiv_py[i]->DrawCopy("E same");
  }
  for(int i=1;i<npriors;i++){
	hdiv[i]->DrawCopy("E same");
  }
  
  
  TLegend *legspectra = new TLegend(0.6, 0.70, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  for(int i=0;i<npriors_py;i++){
  legspectra->AddEntry(hdiv_py[i], Form("pythia w/ p_{T}^{leading} cut (R=%.1lf) / %s",R[i],prior_type_name[priors[0]].Data()), "lp");
  }
  for(int i=1;i<npriors;i++){
  legspectra->AddEntry(hdiv[i], Form("%s/%s",prior_type_name[priors[i]].Data(),prior_type_name[priors[0]].Data()), "lp");
  }
  legspectra->DrawClone("same");
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.4, 0.8,Form("p_{T}^{leading}>%.1lf GeV/c",pTthresh));
    
  str = Form("prior_ratios_pTlead%.0lf_%s.gif",pTthresh,suff.Data());
  cratio->SaveAs(str.Data());
  
  TCanvas *cratio_py = new TCanvas("cratio_py","cratio_py",10,10,can_x,can_y);
  cratio_py->cd();
  cratio_py->SetGrid();
  cratio_py->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T, corr}^{charged} (GeV/c)");
  frame->GetYaxis()->SetTitle("ratio");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(1E-3, 1E4);
  str="Biased PYTHIA over unbiased PYTHIA (particle level)";
  frame->SetTitle(str);
  frame->DrawCopy("");
  for(int i=0;i<npriors_py;i++){
	hdiv_py0[i]->DrawCopy("E same");
  }

    TLegend *legspectra = new TLegend(0.6, 0.70, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  for(int i=0;i<npriors_py;i++){
  legspectra->AddEntry(hdiv_py0[i], Form("PYTHIA: p_{T}^{leading}>%.0lf/p_{T}^{leading}>0 (R=%.1lf)",pTthresh,R[i]),"lp");
  }
  legspectra->DrawClone("same");
  
  str = Form("pythia_prior_ratios_pTlead%.0lf_%s.gif",pTthresh,suff.Data());
  cratio_py->SaveAs(str.Data());
}