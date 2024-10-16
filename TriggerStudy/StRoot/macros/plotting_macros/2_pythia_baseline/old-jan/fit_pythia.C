Double_t TsalisFitFunc(Double_t* x_val, Double_t* par)
{
   double A=par[0];
   double n=par[1];
   double T=par[2];
   double pT=x_val[0];
   
   double y=A*pT*TMath::Power((1.+pT/(n*T)),-n);
   return y;
}


Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0, mu;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
	 mu   = par[4];
    pT=x_val[0];
	 
    Double_t mT = TMath::Sqrt((pT-mu)*(pT-mu)+m0*m0);
	 	 
	y = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
	 return y;
}

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

fit_matt(Double_t R=0.3, Double_t pTcut=5)
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
  
  TString trigger="MB";
  
    Float_t spectraXmin=0;
  Float_t spectraXmax=50;
  Float_t spectraYmin=1E-10;
  Float_t spectraYmax=1E-3;
  

  TString str = "PythiaHistosChargedJets2.root";
  TFile *f = new TFile(str.Data(), "OPEN");


  str=Form("h_PythiaJetPt0p%.0lfTotpTLeading%.0lf",R*10,pTcut);
  if(pTcut<0.1)str=Form("h_PythiaJetPt0p%.0lfTot",R*10);
  TH1D* hPythia=(TH1D*) f->Get(str);
  hPythia->SetLineColor(kBlack);
  hPythia->SetMarkerStyle(29);
  hPythia->SetLineWidth(2);
  hPythia->SetMarkerColor(kBlack);
  hPythia->SetMarkerSize(1.0);
//hPythia->Scale(10);

   TF1* fjet1 = new TF1("fjet1",LevyFitFunc,0.0,50.0,5);
	fjet1->SetParameter(0,1); // B, changes the amplitude, 0.1
	fjet1->SetParameter(1,0.99); // T, changes slope, 0.4
	fjet1->SetParameter(2,18.1); // n, changes how fast spectrum drops, 5.8
	fjet1->SetParameter(3,0.000001); // m0, changes the width, 0.0001		
	fjet1->SetParameter(4,-4.7); // mu, changes the x-axis shift
	fjet1->SetRange(pTcut,50);
	fjet1->SetLineWidth(2);
	fjet1->SetLineStyle(1);
	fjet1->SetLineColor(kGreen);
	fjet1->SetNpx(1000);
  	
	hPythia->Fit("fjet1","R");
  
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
  str="Toymodel hard jet spectrum";
  frame->SetTitle(str);
  frame->DrawCopy("");
  hPythia->DrawCopy("E same");
  fjet1->Draw("same");
  
  
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(hPythia,Form("Pyt. ch. jets w/ p_{T}^{lead}>%.0lfGeV/c, R=%.1lf",pTcut,R),"lp");
  legspectra->AddEntry(fjet1,"fit","l");
  legspectra->DrawClone("same");
  
}

void fit_jan(Float_t R=0.3, Float_t pTcut=5, TString suffix="")
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
  
  TString trigger="MB";
  
    Float_t spectraXmin=0;
  Float_t spectraXmax=50;
  Float_t spectraYmin=1E-10;
  Float_t spectraYmax=1E-3;
  
	Float_t ppxSection=42;
	  
  
  TString str = Form("histos_pythiajet_R%.1lf%s.root",R,suffix.Data());
  TFile *f = new TFile(str.Data(), "OPEN");

  str=Form("hpT_pTl%.0lf",pTcut);
  TH1D* hPythia=(TH1D*) f->Get(str);
  hPythia->SetLineColor(kBlack);
  hPythia->SetMarkerStyle(29);
  hPythia->SetLineWidth(2);
  hPythia->SetMarkerColor(kBlack);
  hPythia->SetMarkerSize(1.0);
  
  TH1I* heventsp=(TH1I*) f->Get("hevts");
	/*Int_t neventsp=heventsp->GetEntries();
   hPythia->Scale(ppxSection/neventsp);*/
   Int_t neventsp=(Int_t) heventsp->GetEntries()/14; //14 pT hard bins
   cout<<"n events:"<<neventsp<<endl;
   hPythia->Scale(1.0/neventsp);
   
   TF1* fjet1 = new TF1("fjet1",PythiaFitFunc,0.0,50.0,10);
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
	/*fjet1->SetParameter(5,0.1);//Amplitude of power law distribution   
   fjet1->SetParameter(6,4); //power*/
	fjet1->SetRange(pTcut,50);
	fjet1->SetLineWidth(2);
	fjet1->SetLineStyle(1);
	fjet1->SetLineColor(kGreen);
	fjet1->SetNpx(1000);
  	
   
   TF1* fjet2 = new TF1("fjet2",TsalisFitFunc,0.0,50.0,3);
   //DEFAULT PARAMETERS: full pythia, R=0.6: [3,10.7,0.44]
   fjet2->SetParameter(0,0.009); // A
   fjet2->SetParameter(1,15); // n
   fjet2->SetParameter(2,0.9); // T
   fjet2->SetRange(pTcut,50);
   fjet2->SetLineWidth(2);
   fjet2->SetLineStyle(1);
   fjet2->SetLineColor(kMagenta);
   fjet2->SetNpx(1000);
   
   
   
   TF1* fjet3 = new TF1("fjet3",LevyFitFunc,0.0,50.0,5);
   fjet3->SetParameter(0,0.1); // B, changes the amplitude, 0.1
   fjet3->SetParameter(1,0.9); // T, changes slope, 0.4
   fjet3->SetParameter(2,5); // n, changes how fast spectrum drops, 5.8
   fjet3->SetParameter(3,0.0001); // m0, changes the width, 0.0001    
   fjet3->SetParameter(4,0); // mu, changes the x-axis shift
   fjet3->SetRange(pTcut,50);
   fjet3->SetLineWidth(2);
   fjet3->SetLineStyle(1);
   fjet3->SetLineColor(kBlue);
   fjet3->SetNpx(1000);
   
	hPythia->Fit("fjet2");
   //hPythia->Fit("fjet3");
   //fjet1->SetParameter(2,10.9);
   //fjet1->SetParameter(7,28.5);
   //fjet1->SetParameter(6,1.99);
   
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
  str="Toymodel hard jet spectrum";
  frame->SetTitle(str);
  frame->DrawCopy("");
  hPythia->DrawCopy("E same");
  fjet2->Draw("same");
  
  
  TLegend *legspectra = new TLegend(0.6467, 0.60, 0.89, 0.90);
  legspectra->SetTextSize(0.03);
  legspectra->SetFillStyle(0);
  legspectra->SetBorderSize(0);
  legspectra->AddEntry(hPythia,Form("Pyt. ch. jets w/ p_{T}^{lead}>%.0lfGeV/c, R=%.1lf",pTcut,R),"lp");
  legspectra->AddEntry(fjet1,"fit","l");
  //legspectra->AddEntry(fjet1,"fit","l");
  legspectra->DrawClone("same");
  
  //cout<<"Integral"<<fjet1->Integral(3,50)<<endl;
}