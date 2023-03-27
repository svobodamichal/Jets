Double_t PythiaFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y,y1,y2, B1, T1, n1, m1, mu1,B2, T2, n2, m2, mu2;
    B1    = par[0];
    T1    = par[1];
    n1    = par[2];
    m1   = par[3];
	 mu1  = par[4];
	 B2    = par[5];
    T2    = par[6];
    n2    = par[7];
    m2   = par[8];
	 mu2   = par[9];

	 pT=x_val[0];
	 
    Double_t mT1 = TMath::Sqrt((pT-mu1)*(pT-mu1)+m1*m1);
	 Double_t mT2 = TMath::Sqrt((pT-mu2)*(pT-mu2)+m2*m2);
	 
	y1 =B1/TMath::Power(1.0+(mT1-m1)/(n1*T1),n1);
	y2=B2/TMath::Power(1.0+(mT2-m2)/(n2*T2),n2);
	
	
	
	if(pT<15)
		y=y1;
	else if(pT<25)
	{
		Double_t c1=(pT-15.0)/10.0;
		y=((1-c1)*y1+c1*y2);
	}
	 else y=y2;

	 return y;
}

Double_t HadronFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B1, T1, n1, m1, mu1;
    B1    = par[0];
    T1    = par[1];
    n1    = par[2];
    m1   = par[3];
	 mu1  = par[4];
	 
	 pT=x_val[0];
	 
    Double_t mT1 = TMath::Sqrt((pT-mu1)*(pT-mu1)+m1*m1);
	 
	y =B1/TMath::Power(1.0+(mT1-m1)/(n1*T1),n1);
	return y;
}

Double_t PhenixHadron(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, A, p0, n;
    A    = par[0];
    p0    = par[1];
    n    = par[2];
    
	 
	 pT=x_val[0];
	 
	 Double_t pTmax=1.6;
	 Double_t R=1.6;
	 Double_t r=(pT>pTmax)? R : R-0.28*(pTmax-pT)*(pTmax-pT);
    Double_t dcross = A*TMath::Power(p0/(p0+pT),n)*r;

	return dcross;
}
/*
void plot_my_plots(TString suffix="_effionly")
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
  
  
   //canvas size
   Int_t can_x=1200; //1600
   Int_t can_y=680; //900
    
  Float_t spectraXmin=0;
  Float_t spectraXmax=50;
  Float_t spectraYmin=1E-12;
  Float_t spectraYmax=1E0;
  
	TString name=Form("histos_pythiajet_R0.3%s.root",suffix.Data());
	TFile *infile = new TFile(name.Data(),"OPEN");
	TH1D* hptl0=(TH1D*) infile->Get("hpT_pTl0");
	TH1D* hptl5=(TH1D*) infile->Get("hpT_pTl5");
	TH1D* hptl5_dete=(TH1D*) infile->Get("hpT_pTl5_dete");
	
	TH1I* hevents=(TH1I*) infile->Get("hevts");
	Int_t nevents=hevents->GetEntries();
	
	name="histos_pythiajet_R0.3_effi_pTsmear.root";
	TFile *infile2 = new TFile(name.Data(),"OPEN");
	TH1D* hptl5_dete2=(TH1D*) infile2->Get("hpT_pTl5_dete");

	
	hptl0->Scale(1./nevents);
	hptl5->Scale(1./nevents);
	hptl5_dete->Scale(1./nevents);
	hptl5_dete2->Scale(1./nevents);

	hptl0->Rebin(4);
	hptl5->Rebin(4);
	hptl5_dete->Rebin(4);
	hptl5_dete2->Rebin(4);
	
	hptl0->SetLineColor(kBlack);
	hptl5->SetLineColor(kRed);
	hptl5_dete->SetLineColor(kBlue);
	hptl5_dete2->SetLineColor(kGreen+1);
	
	hptl0->SetLineWidth(2);
	hptl5->SetLineWidth(2);
	hptl5_dete->SetLineWidth(2);
	hptl5_dete2->SetLineWidth(2);

	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);

   TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
   cspectra->cd();
   cspectra->SetGrid();
   cspectra->SetLogy();
   frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   frame->GetYaxis()->SetTitle("1/N_{events} d^{2}N/dp_{T}d#eta 1/2#pi (1/GeV)");
   frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
   frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
   //frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
   TString str="PYTHIA inclusive jet spectrum";
	//TString str="PYTHIA jet spectrum - highest p_{T} jets only";

   frame->SetTitle(str);
   frame->DrawCopy("");
   hptl0->Draw("same");
   hptl5->Draw("same");
   hptl5_dete->Draw("same");
   hptl5_dete2->Draw("same");

  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(hptl0,"p_{T}^{lead}>0GeV","l");
  legspectra->AddEntry(hptl5,"p_{T}^{lead}>5GeV","l");
  legspectra->AddEntry(hptl5_dete,"p_{T}^{lead}>5GeV, trk. eff. applied","l");
  legspectra->AddEntry(hptl5_dete2,"p_{T}^{lead}>5GeV, trk. eff. and p_{T} smearing","l");
  legspectra->DrawClone("same");
}

void compare_with_matt(Float_t R=0.3, Float_t pTlead=5, TString suffix="")
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
  
  
   //canvas size
   Int_t can_x=1200; //1600
   Int_t can_y=680; //900
    
  Float_t spectraXmin=0;
  Float_t spectraXmax=50;
  Float_t spectraYmin=1E-13;
  Float_t spectraYmax=1E-1;
  
  //my PYTHIA
	TString name=Form("histos_pythiajet_R%.1lf%s.root",R,suffix.Data());
	TFile *infile = new TFile(name.Data(),"OPEN");
	TH1D* hptl=(TH1D*) infile->Get(Form("hpT_pTl%.0lf",pTlead));
	hptl->SetLineColor(kRed);
	hptl->SetMarkerStyle(29);
   hptl->SetLineWidth(2);
   hptl->SetMarkerColor(kRed);
	hptl->SetMarkerSize(0.9);
	
	TH1I* hevents=(TH1I*) infile->Get("hevts");
 	Int_t nevents=hevents->GetEntries()/14.0; //14 pThardbins
	
	Float_t ppXsection=42;
	Float_t jetscale=1.0/nevents;//ppXsection/(nevents);
   //hptl->Rebin(4);
	hptl->Scale(jetscale);
	

  TString str = "PythiaHistosChargedJets2.root";
  TFile *f = new TFile(str.Data(), "OPEN");

//MATT PYTHIA
  str=Form("h_PythiaJetPt0p%.0lfTotpTLeading%.0lf",R*10,pTlead);
  TH1D* hPythia=(TH1D*) f->Get(str);
  hPythia->SetLineColor(kBlack);
  hPythia->SetMarkerStyle(29);
  hPythia->SetLineWidth(2);
  hPythia->SetMarkerColor(kBlack);
  hPythia->SetMarkerSize(1.0);
  
  TH1D* hdiv=(TH1D*) hptl->Clone("hdiv");
  hdiv->Reset("MICE");
  
  
  
  for(int pT=1;pT<50;pT++)
  {
	  
	  int bin1=hPythia->GetXaxis()->FindBin(pT);
	  int bin2=hptl->GetXaxis()->FindBin(pT);
	  if(hPythia->GetBinContent(bin1)==0)continue;
	  float div=hptl->GetBinContent(bin2)/hPythia->GetBinContent(bin1);
	  cout<<"pT :"<<pT<<" jan/matt: "<<div<<endl;
  }
  
  
  
  TH1 *frame = new TH1I("frame", "", 1000, -100, +100);

  TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
  cspectra->cd();
  cspectra->SetGrid();
  cspectra->SetLogy();
  frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  frame->GetYaxis()->SetTitle("d#sigma/dp_{T} (mb/GeV)");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  //frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
  str=Form("PYTHIA charged jets w/ p_{T}^{lead}>%.0lfGeV/c, R=%.1lf",pTlead,R);
  frame->SetTitle(str);
  frame->DrawCopy("");
  hPythia->DrawCopy("E same");
  hptl->Draw("E same");
  
    TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(hPythia,"Matt","lp");
  legspectra->AddEntry(hptl,"Jan","lp");
  legspectra->DrawClone("same");
  
    TString str=Form("obr/PYTHIA_Jan_vs_Matt_R%.1lf_pTlead%0.lf.gif",R,pTlead);
  cspectra->SaveAs(str.Data());
  
}
*/
void compare_with_star(bool star2009=0/*use 2004 or 2009 data*/, Float_t pTlead=0, /*TString suffix="_full",*/TString ext="pdf")
{
	Float_t R=0.4;
	if(star2009)R=0.6;
	
	bool showPytRat=0; //show also ratio of PYTHIA to PYTHIA fit
	
	Color_t cSTAR04=kBlue;
	Color_t cSTAR09=kBlue;
	
	gStyle->SetOptStat(0);
   gStyle->SetPalette(1);
   gStyle->SetOptDate(0);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.09);
   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetTitleSize(0.075,"Y");
   gStyle->SetTitleOffset(0.85,"Y");
   gStyle->SetTitleSize(0.075,"X");
   gStyle->SetTitleOffset(0.95,"X");
   gStyle->SetLabelSize(0.05,"X");
   gStyle->SetLabelSize(0.05,"Y");
  
   //canvas size
   Int_t can_x=1200; //1600
   Int_t can_y=680; //900
    
  Float_t spectraXmin=0;
  Float_t spectraXmax=50;
  Float_t spectraYmin=1E-10;
  Float_t spectraYmax=1E-2;
	Float_t ratioYmin=0.2;
	Float_t ratioYmax=3.0;
	
  Float_t pythiaXmin=7.0;
  
  float ratYmax=1.5;
  float ratYmin=0.5;
  float fitXmin=10;
  float fitXmax=50;
  if(!star2009) 
  {
	  ratYmax=2.5;
	  ratYmin=0.3;
	  fitXmin=5;
  }
  
  TString outdir="./obr";
  
   TF1* fjet1 = new TF1("fjet1",PythiaFitFunc,pythiaXmin,50.0,10);
	fjet1->SetParameter(0,1); // B, changes the amplitude, 0.1
	fjet1->SetParameter(1,0.99); // T, changes slope, 0.4
	fjet1->SetParameter(2,18.1); // n, changes how fast spectrum drops, 5.8
	fjet1->SetParameter(3,0.000001); // m0, changes the width, 0.0001		
	fjet1->SetParameter(4,-4.7); // mu, changes the x-axis shift
	fjet1->SetParameter(5,1); // B, changes the amplitude, 0.1
	fjet1->SetParameter(6,0.99); // T, changes slope, 0.4
	fjet1->SetParameter(7,18.1); // n, changes how fast spectrum drops, 5.8
	fjet1->SetParameter(8,0.000001); // m0, changes the width, 0.0001		
	fjet1->SetParameter(9,-4.7); // mu, changes the x-axis shift
	fjet1->SetRange(fitXmin,fitXmax);
	fjet1->SetLineWidth(2);
	fjet1->SetLineStyle(1);
	fjet1->SetLineColor(kMagenta);
	fjet1->SetNpx(1000);
  
  //my PYTHIA
	//TString name=Form("histos_pythiajet_R%.1lf%s.root",R,suffix.Data());
	TString name=Form("../../plotting_out/root/pp_Fuqiang/Jan_fastjet3/starjet_pythia_R%.1lf_charged0.root",R);
	TFile *infile = new TFile(name.Data(),"OPEN");
	TH1D* hptl=(TH1D*) infile->Get(Form("hpT_pTl%.0lf",pTlead));

	
	TH1I* hevents=(TH1I*) infile->Get("hevts");
	Int_t nevents=hevents->GetEntries()/14.0; //14 pThardbins
	
	Float_t ppXsection=42;
	Float_t jetscale=1.0/nevents; //ppXsection/(nevents);
  	hptl->Scale(jetscale);
	hptl->Rebin(8);
	hptl->Scale(1./8);
	hptl->Fit("fjet1","R");
	
	hptl->SetLineColor(kRed);
	hptl->SetMarkerStyle(22);
   hptl->SetLineWidth(2);
   hptl->SetMarkerColor(kRed);
	hptl->SetMarkerSize(1.0);
	
	//create a white histogram to hide points below pythiaXmin
	TH1D* hptl_cover=hptl->Clone("hptl_cover");
	hptl_cover->Reset("MICE");
	for(int bn=1; bn<hptl_cover->GetNbinsX()+1;bn++)
	{
		double pTc=hptl_cover->GetBinCenter(bn);
		if(pTc<pythiaXmin) 
		{
			hptl_cover->SetBinContent(bn,hptl->GetBinContent(bn)); 
			hptl_cover->SetBinError(bn,hptl->GetBinError(bn));
		}
	}
	hptl_cover->SetLineColor(kWhite);
	hptl_cover->SetMarkerStyle(22);
   hptl_cover->SetLineWidth(3);
   hptl_cover->SetMarkerColor(kWhite);
	hptl_cover->SetMarkerSize(2.0);
	
	//pythia over fit
	TH1D* hpyt_rat=(TH1D*) hptl->Clone("hpyt_rat");
	for(int bn=1; bn<=hpyt_rat->GetNbinsX(); bn++)
	{
		float xpyt=hpyt_rat->GetBinCenter(bn);
		float ypyt=hpyt_rat->GetBinContent(bn);
		
		if (xpyt<fitXmin || xpyt>fitXmax) hpyt_rat->SetBinContent(bn,0);
		else 
		{
			float yfit=fjet1->Eval(xpyt);
			hpyt_rat->SetBinContent(bn,ypyt/yfit);
			hpyt_rat->SetBinError(bn,hpyt_rat->GetBinError(bn)/yfit);
		}
	}
	
	TH1 *frame = new TH1I("frame", "", 1000, -100, +100);
 
	TCanvas *spectra = new TCanvas("cspectra","cspectra",10,10,can_x,can_y);
   cspectra->cd();
   cspectra->SetGrid();
	
   Double_t eps=0.02;
   TPad* p1 = new TPad("p1","p1",0,0.35-eps,1,1,0); p1->Draw();
   p1->SetBottomMargin(eps);
	//p1->SetGrid();
	p1->SetLogy();
	
	TPad* p2 = new TPad("p2","p2",0,0,1,0.35*(1.-eps),0); p2->Draw(); 
	p2->SetTopMargin(0);
	p2->SetBottomMargin(0.35);
	//p2->SetGrid();
	p2->SetFillColor(0);
   p2->SetFillStyle(0);
	
   p1->cd();
	
  //cspectra->SetLogy();
  //frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //frame->GetYaxis()->SetTitle("d#sigma/dp_{T} (#mub/GeV)");
  frame->GetXaxis()->SetRangeUser(spectraXmin, spectraXmax);
  frame->GetYaxis()->SetRangeUser(spectraYmin, spectraYmax);
  //frame->GetYaxis()->SetRangeUser(1E-3, 1E5);
  TString str="";//Form("PYTHIA vs STAR, full jets, R=%.1lf",R);
  //frame->SetTitle(str);
  frame->DrawCopy("");
  fjet1->Draw("same");
  
  //STAR 2009 pp data
  const Int_t nbins=10;
  double avalx[nbins]={10.7032,12.6309,14.9538,17.532,20.5296,24.2735,28.6826,33.7055,39.6978,46.6486};
  double avaly[nbins]={1.81048,0.518655,0.180917,0.0560875,0.0170113,0.00477476,0.00123286,0.000272281, 5.28999e-05,7.37398e-06};

  double aexl[nbins];
  double aexh[nbins];
  double aeyl_stat[nbins]={0.0754247,0.0141737,0.00341262,0.00121029,0.000238343,5.79578e-05,1.5516e-05,3.97909e-06,1.17925e-06,3.11505e-07};
  double aeyh_stat[nbins]={0.0754247,0.0141737,0.00341262,0.00121029,0.000238343,5.79578e-05,1.5516e-05,3.97909e-06,1.17925e-06,3.11505e-07};
  double aeyl_sys[nbins]={0.0876499,0.0402103,0.0175924,0.00620502,0.00214948,0.000574471,0.000188659,4.60589e-05,1.02032e-05,1.57978e-06};
  double aeyh_sys[nbins]={0.0876499,0.0402103,0.0175924,0.00620502,0.00214948,0.000574471,0.000188659,4.60589e-05,1.02032e-05,1.57978e-06};
  
  for(int bin=0; bin<nbins;bin++)
	  {

		  aexl[bin]=0.2;
		  aexh[bin]=0.2;
		  avaly[bin]=avaly[bin]/(1000*2*TMath::Pi());
		  aeyl_stat[bin]=aeyl_stat[bin]/(1000*2*TMath::Pi());
		  aeyh_stat[bin]=aeyh_stat[bin]/(1000*2*TMath::Pi());
		  aeyl_sys[bin]=aeyl_sys[bin]/(1000*2*TMath::Pi());
		  aeyh_sys[bin]=aeyh_sys[bin]/(1000*2*TMath::Pi());
	  }

  
  //STAR 2006 pp data
   const Int_t nbins2=9;
  double bvalx[nbins2]={8.3,10.3,12.6,15.5,19.0,23.4,28.7,35.3,43.3};
  double bvaly[nbins2]={6.4e05,2.4e05,5.9e04,1.19e04,3.5e03,7.8e02,1.41e02,2.2e01,2.6e00};
  double bexl[nbins2]={8.3-7.6,10.3-9.3,12.6-11.4,15.5-14.1,19.0-17.3,23.4-21.3,28.7-26.2,35.3-32.2,43.3-39.6};
  double bexh[nbins2]={9.3-8.3,11.4-10.3,14.1-12.6,17.3-15.5,21.3-19.0,26.2-23.4,32.2-28.7,39.6-35.3,48.7-43.3};
  double beyl_stat[nbins2]={0.9e05,0.3e05,0.6e04,0.1e04,0.2e03,0.4e02,0.09e02,0.3e01,0.6e00};
  double beyh_stat[nbins2]={0.9e05,0.3e05,0.6e04,0.1e04,0.2e03,0.4e02,0.09e02,0.3e01,0.6e00};
  double beyl_sys[nbins2]={3.1e05,1.1e05,2.8e04,0.57e04,1.7e03,3.7e02,0.68e02,1.1e01,1.3e00};
  double beyh_sys[nbins2]={3.1e05,1.1e05,2.8e04,0.57e04,1.7e03,3.7e02,0.68e02,1.1e01,1.3e00};
	
  for(int bin=0; bin<nbins2;bin++)
	  {

		  bvaly[bin]=bvaly[bin]/1E9;
		  beyl_stat[bin]=beyl_stat[bin]/1E9;
		  beyh_stat[bin]=beyh_stat[bin]/1E9;
		  beyl_sys[bin]=beyl_sys[bin]/1E9;
		  beyh_sys[bin]=beyh_sys[bin]/1E9;
	  }
	  
	    
  //STAR pp vs PYTHIA fit
  int nbins3tmp = (star2009)  ? nbins : nbins2;
  const int nbins3 = nbins3tmp;
  double avaly_vs_fit[nbins3];
  double aeyl_stat_vs_fit[nbins3];
  double aeyh_stat_vs_fit[nbins3];
  double aeyl_sys_vs_fit[nbins3];
  double aeyh_sys_vs_fit[nbins3];
  if(star2009)
  for(int bn=0;bn<nbins3;bn++)
  {
	  avaly_vs_fit[bn]=avaly[bn]/fjet1->Eval(avalx[bn]);
	  aeyl_stat_vs_fit[bn]=aeyl_stat[bn]/fjet1->Eval(avalx[bn]);
	  aeyh_stat_vs_fit[bn]=aeyh_stat[bn]/fjet1->Eval(avalx[bn]);
	  aeyl_sys_vs_fit[bn]=aeyl_sys[bn]/fjet1->Eval(avalx[bn]);
	  aeyh_sys_vs_fit[bn]=aeyh_sys[bn]/fjet1->Eval(avalx[bn]);
  }  
  else
  for(int bn=0;bn<nbins3;bn++)
  {
	  avaly_vs_fit[bn]=bvaly[bn]/fjet1->Eval(bvalx[bn]);
	  aeyl_stat_vs_fit[bn]=beyl_stat[bn]/fjet1->Eval(bvalx[bn]);
	  aeyh_stat_vs_fit[bn]=beyh_stat[bn]/fjet1->Eval(bvalx[bn]);
	  aeyl_sys_vs_fit[bn]=beyl_sys[bn]/fjet1->Eval(bvalx[bn]);
	  aeyh_sys_vs_fit[bn]=beyh_sys[bn]/fjet1->Eval(bvalx[bn]);
  }  
	
	TGraphAsymmErrors* gppjets06 = new TGraphAsymmErrors(nbins2, bvalx, bvaly, bexl, bexh, beyl_stat, beyh_stat);
   gppjets06->SetLineColor(cSTAR04);
   gppjets06->SetLineWidth(2);
	gppjets06->SetMarkerColor(cSTAR04);
	gppjets06->SetMarkerStyle(33);
	gppjets06->SetMarkerSize(2.0);
	gppjets06->GetXaxis()->SetLimits(spectraXmin,spectraXmax);
	gppjets06->GetHistogram()->SetMaximum(spectraYmax);   // along          
   gppjets06->GetHistogram()->SetMinimum(spectraYmin);
	//gppjets06->SetTitle(Form("R_{AA}, R=%.1lf, p_{T}^{leading}>%.1lf GeV/c",R,pTthresh));
	gppjets06->SetTitle(str);
	gppjets06->GetXaxis()->SetTitle("p_{T,jet}^{full} (GeV/c)");
   gppjets06->GetYaxis()->SetTitle("d^{2}#sigma/(dp_{T}d#eta2#pi) (mb/GeV) ");
	//gppjets06->GetXaxis()->CenterTitle()
	gppjets06->GetYaxis()->CenterTitle();
	gppjets06->GetYaxis()->SetTitleOffset(0.8);

	
	TGraphAsymmErrors* gppjet06_sys = new TGraphAsymmErrors(nbins2, bvalx, bvaly, bexl, bexh, beyl_sys, beyh_sys);
	gppjet06_sys->SetFillColor(cSTAR04);
   gppjet06_sys->SetLineColor(cSTAR04);
	gppjet06_sys->SetLineWidth(1);
   gppjet06_sys->SetFillStyle(0);

	if(!star2009) 
	{
		gppjets06->Draw("ap");
		gppjet06_sys->Draw("2");
	}
	

	  
  TGraphAsymmErrors* gppjets09 = new TGraphAsymmErrors(nbins, avalx, avaly, aexl, aexh, aeyl_stat, aeyh_stat);
   gppjets09->SetLineColor(cSTAR09);
   gppjets09->SetLineWidth(2);
	gppjets09->SetMarkerColor(cSTAR09);
	gppjets09->SetMarkerStyle(29);
	gppjets09->SetMarkerSize(2.0);
	gppjets09->GetXaxis()->SetLimits(spectraXmin,spectraXmax);
	gppjets09->GetHistogram()->SetMaximum(spectraYmax);   // along          
   gppjets09->GetHistogram()->SetMinimum(spectraYmin);
	//gppjets09->SetTitle(Form("R_{AA}, R=%.1lf, p_{T}^{leading}>%.1lf GeV/c",R,pTthresh));
	gppjets09->SetTitle(str);
	gppjets09->GetXaxis()->SetTitle("p_{T,jet}^{full} (GeV/c)");
   gppjets09->GetYaxis()->SetTitle("d^{2}#sigma/(dp_{T}d#eta2#pi) (mb/GeV) ");
	//gppjets09->GetXaxis()->CenterTitle()
	gppjets09->GetYaxis()->CenterTitle();
	gppjets09->GetYaxis()->SetTitleOffset(0.8);

	
	TGraphAsymmErrors* gppjet09_sys = new TGraphAsymmErrors(nbins, avalx, avaly, aexl, aexh, aeyl_sys, aeyh_sys);
	gppjet09_sys->SetFillColor(cSTAR09);
   gppjet09_sys->SetLineColor(cSTAR09);
	gppjet09_sys->SetLineWidth(1);
   gppjet09_sys->SetFillStyle(0);
	
	if(star2009) 
	{
		gppjets09->Draw("ap");
		gppjet09_sys->Draw("2");
	}
	  
	  
  	hptl->Draw("histo E same");
	hptl_cover->Draw("histo E same"); //hide points below pythiaXmin
	
	TLine *one = new TLine(spectraXmin, spectraYmax, pythiaXmin, spectraYmax);
	one->SetLineWidth(1);
	one->SetLineStyle(1);
	one->SetLineColor(kBlack);
	one->DrawClone("same");

	//plot the STAR data once again so they are above PYTHIA
	if(star2009)
		gppjets09->Draw("p");
	else
		gppjets06->Draw("p");


	 double posx1=0.48;
    double posx2=0.90;
    double posy1=0.65;
    double posy2=0.85;
	
	TLegend *legraa = new TLegend(posx1,posy1,posx2,posy2);
	legraa->SetTextSize(0.065);
	legraa->SetFillStyle(0);
	legraa->SetBorderSize(0);
	if(star2009) 
		legraa->AddEntry(gppjets09, "STAR, anti-k_{T}","lp");
	else
		legraa->AddEntry(gppjets06, Form("STAR, midpoint-cone, R=%.1lf",R), "lp");
	legraa->AddEntry(hptl, Form("PYTHIA 6, anti-k_{T}, R=%.1lf",R), "lp");
   // legraa->AddEntry(gtaa, "T_{AA} uncertainty", "f");
	legraa->DrawClone("same");
	
	TLegend *model_info = new TLegend(0.2, 0.15, 0.35, 0.45);
	model_info->SetTextSize(0.065);
	model_info->SetFillStyle(0);
	model_info->SetBorderSize(0);
	model_info->SetMargin(0.05);
	model_info->SetHeader("p+p #sqrt{s_{NN}}=200 GeV");
	//model_info->AddEntry("",Form("R=%.1lf",R),"");
	model_info->DrawClone("same");
	
	//RATIO plot
	//*****************************
	p2->cd();
	
	frame->GetYaxis()->SetRangeUser(ratioYmin, ratioYmax);
	frame->GetYaxis()->SetTitleSize(0.1);
	frame->GetYaxis()->SetTitleOffset(0.5);
   frame->GetYaxis()->SetLabelSize(0.08);
   frame->GetXaxis()->SetTitleSize(0.15);
   frame->GetXaxis()->SetTitleOffset(0.9);
   frame->GetXaxis()->SetLabelSize(0.1);
	frame->GetXaxis()->SetTitle("p_{T,jet}^{full} (GeV/c)");
   frame->GetYaxis()->SetTitle("data/PYTHIA");
	frame->DrawCopy();
	
   TGraphAsymmErrors* gpp_vs_fit;
	if(star2009) gpp_vs_fit = new TGraphAsymmErrors(nbins3, avalx, avaly_vs_fit, aexl, aexh, aeyl_stat_vs_fit, aeyh_stat_vs_fit);
	else gpp_vs_fit = new TGraphAsymmErrors(nbins3, bvalx, avaly_vs_fit, bexl, bexh, aeyl_stat_vs_fit, aeyh_stat_vs_fit);
   gpp_vs_fit->SetLineColor(cSTAR09);
	gpp_vs_fit->SetMarkerColor(cSTAR09);
	gpp_vs_fit->SetMarkerStyle(29);
	if(!star2009)
	{
		gpp_vs_fit->SetLineColor(cSTAR04);
		gpp_vs_fit->SetMarkerColor(cSTAR04);
		gpp_vs_fit->SetMarkerStyle(33);
	}
	
   gpp_vs_fit->SetLineWidth(2);
	gpp_vs_fit->SetMarkerSize(2.0);
	gpp_vs_fit->GetXaxis()->SetLimits(spectraXmin,spectraXmax);
	gpp_vs_fit->GetHistogram()->SetMaximum(ratYmax);   // along          
   gpp_vs_fit->GetHistogram()->SetMinimum(ratYmin);
   gpp_vs_fit->Draw("p");
	//gpp_vs_fit->SetTitle(Form("R_{AA}, R=%.1lf, p_{T}^{leading}>%.1lf GeV/c",R,pTthresh));
	gpp_vs_fit->SetTitle("");
	//gpp_vs_fit->GetXaxis()->CenterTitle()
	gpp_vs_fit->GetYaxis()->CenterTitle();
  gpp_vs_fit->GetYaxis()->SetLabelSize(0.06);
  gpp_vs_fit->GetXaxis()->SetLabelSize(0.06);
  gpp_vs_fit->GetXaxis()->SetTitleSize(0.1);
  gpp_vs_fit->GetYaxis()->SetTitleOffset(0.5);
  gpp_vs_fit->GetYaxis()->SetTitleSize(0.08);
	
	if(star2009) TGraphAsymmErrors* gpp_vs_fit_sys = new TGraphAsymmErrors(nbins3, avalx, avaly_vs_fit, aexl, aexh, aeyl_sys_vs_fit, aeyh_sys_vs_fit);
	else TGraphAsymmErrors* gpp_vs_fit_sys = new TGraphAsymmErrors(nbins3, bvalx, avaly_vs_fit, bexl, bexh, aeyl_sys_vs_fit, aeyh_sys_vs_fit);
   gpp_vs_fit_sys->SetLineColor(cSTAR09);
	if(!star2009)
	{
		gpp_vs_fit_sys->SetLineColor(cSTAR04);
	}
	gpp_vs_fit_sys->SetLineWidth(1);
   gpp_vs_fit_sys->SetFillStyle(0);
	gpp_vs_fit_sys->Draw("2");
	
		TLine *one = new TLine(spectraXmin, 1, spectraXmax, 1);
  one->SetLineWidth(3);
  one->SetLineStyle(2);
  one->SetLineColor(kBlack);
  one->DrawClone("same");
  
  if(showPytRat) hpyt_rat->DrawCopy("same");
  
  str = Form("%s/pythia6_R0%.0lf_vs_star.%s", outdir.Data(),R*10,ext.Data());
  spectra->SaveAs(str.Data());
}


void compare_charged_hadrons(TString ext="png")
{
	//settings
	bool showPhenixhPi0=1; //show phenix pi0 spectrum
	bool showSTARChH=1; //show STAR ch. h. spectrum
	bool showPythiaPi0=1; //show PYTHIA pi0 spectrum
	
	//constants
	const double h_pi0_rat=1.6; //hadron/pi0 ratio in p+p collisions
	const double figx=800;
	const double figy=400;
	const double ppXsection=42;
	
	//colors and markers 
	Color_t cstar=kBlue; //STAR data
	Color_t cphx=kAzure-8; //PHENIX
	Color_t cpythia=kRed-9; //PYTHIA
	int markstar=29;
	int markphx=24;
	int markpythia=20;
	int msizestar=2;
	int msizephx=2;
	int msizepythia=1;

		
	//line at unity
	Color_t lineOneColor=kGray;
	int lineOneWidth=1;
	int lineOneStyle=2;
	
	//graph dimensions
	double xmin=0;
	double xmax=15;
	double ymin=1E-7;
	double ymax=1E1;

	//fit functions
	double fitXmin=0.2;
	double fitXmax=15.0;
	/*TF1* ffit = new TF1("ffit",HadronFitFunc,fitXmin,fitXmax,5);
	ffit->SetParameter(0,1); // B, changes the amplitude, 0.1
	ffit->SetParameter(1,0.4); // T, changes slope, 0.4
	ffit->SetParameter(2,3.0); // n, changes how fast spectrum drops, 5.8
	ffit->SetParameter(3,0.0001); // m0, changes the width, 0.0001		
	ffit->SetParameter(4,1.0); // mu, changes the x-axis shift*/
	TF1* ffit=new TF1("ffit",PhenixHadron, fitXmin, fitXmax,3);
	ffit->SetParameter(0,386); //A [mb/(GeV/c^2)]
	ffit->SetParameter(1,1.219); //p0 [GeV/c]
	ffit->SetParameter(2,9.99); //n
	ffit->SetRange(fitXmin,fitXmax);
	ffit->SetLineWidth(2);
	ffit->SetLineStyle(1);
	ffit->SetLineColor(kMagenta);
	ffit->SetNpx(1000);
	
	TF1* fPHXhadron=new TF1("fPHXhadron",PhenixHadron, fitXmin, fitXmax,3);
	fPHXhadron->SetParameter(0,386); //A [mb/(GeV/c^2)]
	fPHXhadron->SetParameter(1,1.219); //p0 [GeV/c]
	fPHXhadron->SetParameter(2,9.99); //n
	fPHXhadron->SetRange(fitXmin,fitXmax);
	fPHXhadron->SetLineWidth(2);
	fPHXhadron->SetLineStyle(1);
	fPHXhadron->SetLineColor(kBlue);
	fPHXhadron->SetNpx(1000);
	
	//*********************
	//input data
	//*********************
	
	//pi0 spectrum in p+p from Phenix*1.6 -> equivalent to ch. hadron spectrum
	const int nphxpi0=17;
	double phxpi0_x[nphxpi0]={1.215, 1.719, 2.223, 2.726, 3.228, 3.73, 4.232, 4.733, 5.234, 5.735, 6.236, 6.737, 7.452, 8.457, 9.46, 10.861, 13.25};
	double phxpi0_x_err[nphxpi0];
	double phxpi0_y[nphxpi0]={0.3733, 0.06052, 0.01221, 0.003308, 0.0009978, 0.0003385, 0.0001187, 4.726E-05, 2.206E-05, 1.113E-05, 4.999E-06, 3.003E-06, 1.08E-06, 4.853E-07, 1.643E-07, 5.227E-08, 1.19E-08};
	for (int i=0; i<nphxpi0; i++){
		phxpi0_y[i]=phxpi0_y[i]*h_pi0_rat;
		phxpi0_x_err[i]=0.1;
	}
	double phxpi0_y_stat[nphxpi0]={0.006117, 0.001075, 0.0003011, 0.0001181, 5.652E-05, 2.456E-05, 2.906E-06, 1.99E-06, 1.096E-06, 5.005E-07, 3.172E-07, 2.304E-07, 9.53E-08, 5.832E-08, 3.166E-08, 1.169E-08, 4.91E-09};
	double phxpi0_y_sys[nphxpi0]={0.02731, 0.004284, 0.0008644, 0.0002375, 7.278E-05, 2.591E-05, 9.888E-06, 4.001E-06, 1.929E-06, 1.028E-06, 4.759E-07, 2.95E-07, 1.092E-07, 5.253E-08, 1.813E-08, 6.137E-09, 1.894E-09};
	
	TGraphErrors* gphxpi0=new TGraphErrors(17,phxpi0_x,phxpi0_y,0,phxpi0_y_stat);
	TGraphErrors* gphxpi0_sys=new TGraphErrors(17,phxpi0_x,phxpi0_y,phxpi0_x_err,phxpi0_y_sys);
	
	gphxpi0->SetLineColor(cphx);
   gphxpi0->SetLineWidth(1);
	gphxpi0->SetMarkerStyle(markphx);
	gphxpi0->SetMarkerSize(msizephx);
	gphxpi0->SetMarkerColor(cphx);
	gphxpi0_sys->SetLineColor(cphx);
   gphxpi0_sys->SetLineWidth(1);
	gphxpi0_sys->SetFillStyle(0);

	//STAR pp ch. hadron spectrum
	const int nstrchh=32;
	double starChH_x[nstrchh]={0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.49, 2.7, 2.9, 3.16, 3.55, 4.06, 4.7, 5.48, 6.42, 7.43, 8.43, 9.44};
	double starChH_x_low[nstrchh]={0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.09, 0.1, 0.1, 0.16, 0.2, 0.26, 0.3, 0.38, 0.42, 0.43, 0.43, 0.44};
	double starChH_x_hi[nstrchh]={0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.11, 0.1, 0.1, 0.19, 0.25, 0.34, 0.4, 0.52, 0.58, 0.57, 0.57, 0.56};
	double starChH_y[nstrchh]={0.58, 0.35, 0.21, 0.13, 0.079, 0.051, 0.033, 0.022, 0.015, 0.01, 0.0072, 0.0051, 0.0035, 0.0025, 0.0018, 0.0013, 0.00094, 0.0007, 0.00051, 0.00038, 0.00024, 0.00015, 9.6E-05, 5.4E-05, 2.2E-05, 7.8E-06, 2.4E-06, 7.2E-07, 1.9E-07, 7.7E-08, 2.4E-08, 1.2E-08};
	double starChH_y_err[nstrchh]={0.04, 0.02, 0.01, 0.01, 0.004, 0.003, 0.002, 0.001, 0.001, 0.001, 0.0004, 0.0003, 0.0003, 0.0002, 0.0001, 0.0001, 9E-05, 7E-05, 6E-05, 4E-05, 3E-05, 2E-05, 1E-05, 4E-06, 2E-06, 6E-07, 2E-07, 7E-08, 2E-08, 1.2E-08, 6E-09, 4E-09};
	double starChH_y_rat[nstrchh];
	double starChH_y_err_rat[nstrchh];

		for (int i=0; i<nstrchh; i++){
		starChH_y[i]=starChH_y[i]*ppXsection;
		//ratio STAR/PHENIX
		starChH_y_rat[i]=starChH_y[i]/fPHXhadron->Eval(starChH_x[i]);
		starChH_y_err_rat[i]=starChH_y_err[i]/fPHXhadron->Eval(starChH_x[i]);
		}
	
	TGraphAsymmErrors* gstarChH=new TGraphAsymmErrors(nstrchh,starChH_x,starChH_y,starChH_x_low,starChH_x_hi,starChH_y_err,starChH_y_err);
	gstarChH->SetLineColor(cstar);
   gstarChH->SetLineWidth(1);
	gstarChH->SetMarkerStyle(markstar);
	gstarChH->SetMarkerSize(msizestar);
	gstarChH->SetMarkerColor(cstar);
	
	//ratio STAR/Phenix
	TGraphAsymmErrors* gstarChH_rat=new TGraphAsymmErrors(nstrchh,starChH_x,starChH_y_rat,starChH_x_low,starChH_x_hi,starChH_y_err_rat,starChH_y_err_rat);
	gstarChH_rat->SetLineColor(cstar);
   gstarChH_rat->SetLineWidth(1);
	gstarChH_rat->SetMarkerStyle(markstar);
	gstarChH_rat->SetMarkerSize(msizestar);
	gstarChH_rat->SetMarkerColor(cstar);
	
	TString name="../../plotting_out/root/pp_Fuqiang/Jan_fastjet3/starjet_pythia_R0.4_charged1_with_hadron_spectrum.root";
	TFile *infile = new TFile(name.Data(),"OPEN");
	TH1D* hpthadron=(TH1D*) infile->Get("hpt_particle");

	TH1I* hevents=(TH1I*) infile->Get("hevts");
	Int_t nevents=hevents->GetEntries()/14.0; //14 pThardbins
	
	Float_t scale=1.0/nevents; //ppXsection/(nevents);
  	hpthadron->Scale(scale);
	//hpthadron->Rebin(8);
	//hpthadron->Scale(1./8);
	//hpthadron->Fit("ffit","R");
	hpthadron->SetLineColor(cpythia);
	hpthadron->SetMarkerStyle(markpythia);
	hpthadron->SetMarkerSize(msizepythia);
	hpthadron->SetMarkerColor(cpythia);
	
	name="../../plotting_out/root/pp_Fuqiang/Jan_fastjet3/starjet_pythia_R0.2_charged0_pi0spec_and_nojets.root";
	TFile *infile2 = new TFile(name.Data(),"OPEN");
	TH1D* hptpi0=(TH1D*) infile2->Get("hpt_pi0");
	hptpi0->Sumw2();
	TH1I* hevents2=(TH1I*) infile2->Get("hevts");
	Int_t nevents2=hevents2->GetEntries()/14.0; //14 pThardbins
	double scale2=1.0/nevents2; //ppXsection/(nevents);
  	hptpi0->Scale(scale2);
	hptpi0->SetLineColor(cpythia);
	hptpi0->SetMarkerStyle(markpythia+1);
	hptpi0->SetMarkerSize(msizepythia);
	hptpi0->SetMarkerColor(cpythia+1);
		
	TCanvas *cfig=new TCanvas("cfig","comparison",10,10,figx,2*figy);
	cfig->Divide(1,2);
	cfig->cd(1);
	gPad->SetMargin(0.20,0.03,0.20,0.03);
	gPad->SetLogy();
	
	TH1D* htmp1=new TH1D("htmp1","",100,xmin,xmax);
	htmp1->SetXTitle("p_{T} [GeV/c]");htmp1->SetYTitle("#frac{1}{p_{T}}#frac{d#sigma}{dydp_{T}}");
	htmp1->SetTitleOffset(1.2,"x");htmp1->SetTitleOffset(1.1,"y");
	htmp1->SetTitleSize(0.075,"x");htmp1->SetTitleSize(0.08,"y");
	htmp1->SetLabelSize(0.055,"x");htmp1->SetLabelSize(0.055,"y");
	//htmp1->SetNdivisions(505,"x");htmp1->SetNdivisions(505,"y");
	htmp1->SetMinimum(ymin);htmp1->SetMaximum(ymax);
	htmp1->DrawCopy();
	
	TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
	legspectra->SetTextSize(0.045);
	legspectra->SetFillStyle(0);
	legspectra->SetBorderSize(0);
    
	hpthadron->Draw("same");
   legspectra->AddEntry(hpthadron,"Pythia 6","p");

	if(showPythiaPi0)
	{
		hptpi0->DrawCopy("histo same");
		legspectra->AddEntry(hptpi0,"Pythia #pi^{0}","l");
	}
	
	if(showSTARChH) {
		gstarChH->Draw("P");
		legspectra->AddEntry(gstarChH,"STAR h^{+}+h^{-}","p");
	}
	
	if(showPhenixhPi0) {
		gphxpi0_sys->Draw("2");
		gphxpi0->Draw("P");
		legspectra->AddEntry(gphxpi0,"PHENIX #pi^{0}*1.6","p");
		fPHXhadron->DrawCopy("same");
		legspectra->AddEntry(fPHXhadron,"PHENIX fit","l");
	}
	
	legspectra->DrawClone("same");
	
	//Claculate ratios PYTHIA/Phenix 
	TH1D* hratPytPhx=(TH1D*) hpthadron->Clone("hratPytPhx");
	for (int bn=1; bn<=hratPytPhx->GetNbinsX(); bn++)
	{
		double pTc=hratPytPhx->GetBinCenter(bn);
		double y=hratPytPhx->GetBinContent(bn);
		double yphx=fPHXhadron->Eval(pTc);
		double yerr=hratPytPhx->GetBinError(bn);
		hratPytPhx->SetBinContent(bn,y/yphx);
		hratPytPhx->SetBinError(bn,yerr/yphx);
	}
	TH1D* hratPytPi0Phx=(TH1D*) hptpi0->Clone("hratPytPi0Phx");
	for (int bn=1; bn<=hratPytPi0Phx->GetNbinsX(); bn++)
	{
		double pTc=hratPytPi0Phx->GetBinCenter(bn);
		double y=hratPytPi0Phx->GetBinContent(bn);
		double yphx=fPHXhadron->Eval(pTc)/h_pi0_rat; //Phenix Pi0 spectrum
		double yerr=hratPytPi0Phx->GetBinError(bn);
		hratPytPi0Phx->SetBinContent(bn,y/yphx);
		hratPytPi0Phx->SetBinError(bn,yerr/(yphx*pTc));
	}
	
	cfig->cd(2);
	gPad->SetMargin(0.20,0.03,0.20,0.03);
	TH1D* htmp2=new TH1D("htmp2","",100,xmin,xmax);
	htmp2->SetXTitle("p_{T} [GeV/c]");htmp2->SetYTitle("ratio");
	htmp2->SetTitleOffset(1.2,"x");htmp2->SetTitleOffset(1.1,"y");
	htmp2->SetTitleSize(0.075,"x");htmp2->SetTitleSize(0.08,"y");
	htmp2->SetLabelSize(0.055,"x");htmp2->SetLabelSize(0.055,"y");
	//htmp1->SetNdivisions(505,"x");htmp2->SetNdivisions(505,"y");
	htmp2->SetMinimum(0.2);htmp2->SetMaximum(3.6);
	htmp2->DrawCopy();
	
	TLegend *legratio = new TLegend(0.22, 0.65, 0.65, 0.90);
	legratio->SetTextSize(0.045);
	legratio->SetFillStyle(0);
	legratio->SetBorderSize(0);
	legratio->AddEntry(hratPytPhx,"PYTHIA/PHENIX","lp");
	
	hratPytPhx->DrawCopy("same");
	
	if(showPythiaPi0)
	{
		hratPytPi0Phx->DrawCopy("same");
		legratio->AddEntry(hratPytPi0Phx,"PYTHIA #pi^{0}/PHENIX #pi^{0}","lp");
	}
	
	if(showSTARChH){
		gstarChH_rat->Draw("P");
		legratio->AddEntry(gstarChH_rat,"STAR/PHENIX","lp");
	}
	legratio->DrawClone("same");

		
	TLine *one = new TLine(xmin, 1, xmax, 1);
	one->SetLineWidth(2);
	one->SetLineStyle(3);
	one->SetLineColor(kBlack);
	one->DrawClone("same");
	
	cfig->SaveAs(Form("obr/PYTHIA_ch_hadrons_vs_STAR_and_Phenix.%s",ext.Data()));
}
